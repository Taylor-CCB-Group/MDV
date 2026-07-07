# 02 — Flexible table handling in `convert-spatial`

> How SpatialData `tables` become MDV datasources today, why the current single-concat merge
> makes a mess for heterogeneous tables, and the design space for fixing it.

## The problem in one sentence

The converter iterates **every table of every store into one flat list** and does **one
`anndata.concat`**, always producing **exactly one obs datasource + one var datasource** — so a
`grid`/bins table and a `cells` table get row-stacked into a single datasource with a zero-filled
union gene matrix, which is meaningless.

## Current pipeline (grounded)

CLI: `convert-spatial` at [cli.py:172](../../../python/mdvtools/cli.py), handler
`convert_spatial` at [cli.py:192](../../../python/mdvtools/cli.py) → builds
`SpatialDataConversionArgs` ([conversion.py:23](../../../python/mdvtools/spatial/conversion.py))
→ `convert_spatialdata_to_mdv` ([conversion.py:968](../../../python/mdvtools/spatial/conversion.py)).

> Two papercuts noticed in passing: the click positional args are declared
> `output_folder` then `spatialdata_path` but the handler signature is
> `(spatialdata_path, output_folder, …)` — click binds by name, so on the CLI you type
> `convert-spatial OUTPUT_FOLDER SPATIALDATA_PATH`. And `--point-transform` is exposed only on the
> argparse / report entrypoints, not the click CLI.

The merge, `_concat_spatial_tables`
([conversion.py:885](../../../python/mdvtools/spatial/conversion.py)), called at
[conversion.py:1117](../../../python/mdvtools/spatial/conversion.py):

```python
def _concat_spatial_tables(adata_objects):
    return ad_concat(adata_objects, index_unique="_", join="outer", fill_value=0)
```

As pseudo-code:

```
adata_objects = []
for store in discovered_stores:
    for (table_name, adata) in read_zarr(store).tables.items():   # EVERY table, no grouping
        resolve_regions_for_table(adata)         # writes obs.x/y/spatial_region, uns.mdv
        adata.obs["spatialdata_path"] = store_name
        adata.obs["table_name"]       = table_name   # the ONLY per-table provenance
        adata_objects.append(adata)                  # flat list, all shapes mixed
merged = anndata.concat(adata_objects, index_unique="_", join="outer", fill_value=0)
convert_scanpy_to_mdv(output, merged, ...)           # -> ONE obs ds + ONE var ds
```

`convert_scanpy_to_mdv` ([conversions.py:118](../../../python/mdvtools/conversions.py)) then maps
`obs`→`cells`, `var`→`genes`, links them with `add_rows_as_columns_link`, and stores `X` as a
`rows_as_columns` subgroup.

### Element → MDV construct (today)

| SpatialData element | MDV construct |
|---|---|
| `tables` obs | rows of the **single** `cells` datasource |
| `tables` var | rows of the **single** `genes` datasource |
| `tables` X / layers | `rows_as_columns` subgroups (gene scores) on `cells` |
| `tables` obsm (spatial, umap) | `x`/`y` + dim-reduction columns |
| `shapes` (annotated region) | GeoJSON `images/<region>.geo.json` + region entry; coords → obs x/y |
| `labels` | coordinate transform / region only; **no datasource** |
| `points` | **not converted** |
| `images` | `viv_image` per region |
| region / region_key / instance_key | consumed to compute regions; only `spatial_region` text column survives |

## Why it makes a mess

SpatialData's data model (`TableModel`, `get_table_keys` →
`(region, region_key, instance_key)`) says **each table annotates a distinct element / entity
type**. The merge ignores this:

- **Rows are stacked across entity types.** `n_obs_merged = Σ n_obs`. A grid table (raster
  lattice, e.g. Visium HD 2µm bins) and a cells table (segmented cells) end up as one datasource
  whose rows are a mix of grid squares and cells. Row count is meaningless; per-cell analysis is
  polluted by grid rows and vice-versa.
- **The gene axis is force-unioned.** `join="outer"` + `fill_value=0` across unrelated feature
  spaces (e.g. a protein panel vs a transcriptome) produces a huge block that is `0` everywhere a
  feature doesn't belong — statistically and visually misleading. (Outer join is deliberate:
  inner join would collapse disjoint gene sets to zero, comment at conversion.py:887.)
- **Row→element is only a text column.** After merge there is one obs datasource but multiple
  regions; a row relates to its element solely via the soft `obs["spatial_region"]` string, not a
  structural link.
- **Per-table embeddings leak into shared columns.** `compute_x_umap` UMAP/leiden are explicitly
  *not* comparable across tables, yet land in shared `X_umap`/`leiden` columns; leiden categories
  are per-table-prefixed as a workaround (`_prefix_table_leiden_categories`, conversion.py:938).

## Design space

The lever: replace the single flat list + single `_concat_spatial_tables`
([conversion.py:1036–1131](../../../python/mdvtools/spatial/conversion.py)) with a **grouping
step** — `group_tables(adata_objects) -> {group_name: [adata]}` — then loop groups, producing one
`(obs, var)` datasource pair per group via `convert_scanpy_to_mdv` (which already accepts custom
datasource names). A new field on `SpatialDataConversionArgs` selects the policy.

| Option | How | Pros | Cons | Effort |
|---|---|---|---|---|
| **1. One datasource per table** | each `(store, table)` → own obs/var pair, names from table name (`_allocate_table_prefix` exists) | never mixes entities; simplest correct default; native gene sets, no zero-fill; embeddings stay per-table | many datasources for big batches; cross-table compare needs explicit links; default-view template assumes single `cells`/`genes` | Moderate |
| **2. Merge by region / annotated element** | group by `get_table_keys` `region` (or `instance_key` space); concat within group | semantically principled; `cells` across stores merge, `cells`+`grid` stay separate; reuses region resolution | needs stable "same element type" identity across stores; no-attrs tables need a fallback bucket | Moderate |
| **3. Merge by compatible shape (auto)** | concat only tables with same `instance_key` semantics and/or var overlap above a threshold + same kind | automatic; minimizes datasource count while avoiding nonsense | fuzzy heuristics can misclassify; hard to predict; needs override | Higher |
| **4. User-specified grouping config** | CLI/JSON mapping tables → target datasources + per-group options | fully explicit, reproducible, arbitrary edge cases | user burden; verbose for batch; needs discovery tooling | Moderate |

**Recommended shape:** default to **(2) merge-by-region** (or **(1) one-per-table** as the safest
default), with **(4) config override**, and **(3) as an opt-in `--auto-group`**. All four converge
on the same refactor.

### Concrete refactor touch-points

- Merge block → per-group concat + per-group `convert_scanpy_to_mdv`:
  [conversion.py:1099–1131](../../../python/mdvtools/spatial/conversion.py).
- Table-collection loop (attach grouping key while iterating):
  [conversion.py:1057–1077](../../../python/mdvtools/spatial/conversion.py).
- Regions-metadata write (currently single-obs-ds hardcoded — must target the right obs ds per
  group): [conversion.py:1191–1203](../../../python/mdvtools/spatial/conversion.py) +
  `project_helpers.build_spatial_regions_metadata`
  ([project_helpers.py:11](../../../python/mdvtools/spatial/project_helpers.py)).
- New arg on `SpatialDataConversionArgs` + expose on click CLI + argparse.
- Default-view template assumes one obs/var pair: `spatial_view_template.json`,
  `set_default_spatial_image_view` ([project_helpers.py:29](../../../python/mdvtools/spatial/project_helpers.py))
  — needs per-group views or a chosen "primary" datasource.
- Reusable helpers already present: `_allocate_table_prefix` (conversion.py:911),
  `_sanitize_table_prefix` (conversion.py:906), `get_table_keys` (region/instance_key).

### Prior art to mine

- **`annotations.py`** (`patch-spatial-annotations`) already uses `instance_key` as a join key
  ([annotations.py:96](../../../python/mdvtools/spatial/annotations.py)) — the pattern a
  merge-by-region mode would reuse for structural row→element association. (It's "V1": requires
  exactly one table and row-order match.)
- **`xenium.py`** (`convert_xenium_to_mdv`) handles a single Xenium store's elements directly —
  self-described as needing review, but mineable for per-element handling ideas.
- **The user's other conversion scripts** (mentioned as living in other repos — e.g. under
  `~/code/spatialdata`, `~/code/py/sd_sketch`) should be pulled in here for comparison before
  fixing the grouping policy; they likely already encode grouping decisions worth adopting.
- `conversion_report.py` is a batch harness (`python -m mdvtools.spatial.conversion` over a
  directory) — use it to regression-test any new merge modes against the existing corpus under
  `data/spatialdata_conversion_report_*`.

## Relationship to the JS DataLoader (theme 1)

If theme 1 lands, a chunk of this converter's job — flattening tables into MDV columns ahead of
time — could instead happen **at read time in the browser**. The grouping *policy* (which tables
are which entity type, how they link to elements) is still needed either way; it just moves from
"bake into h5" to "describe in datasource metadata the JS loader consumes." Designing the grouping
model cleanly now pays off in both places. See [README.md](README.md) → "the radical vision."
