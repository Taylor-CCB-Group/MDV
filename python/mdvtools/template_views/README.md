# Template views

Infer experimental-design fields from an existing MDV project and write a small set of template views. There is no LLM: roles come from `datasources.json` column names and links.

## What it writes

| View | Contents |
| --- | --- |
| `1 Study overview` | Cohort text, design bars, QC histograms, filters |
| `2 Cell atlas` | UMAP(s) coloured by cell type, tissue, disease |
| `3 Marker genes and signatures` | Feature UMAPs, top-20 marker table, factor-varying dots |
| `4 Tissue and disease context` | Condition stacked bars and abundance |
| `5 RNA-protein concordance` | Side-by-side RNA / protein UMAPs when a protein link exists |

A view is skipped when the project lacks the fields it needs (for example no UMAP, no RNA link, no protein). Existing views are left in place; these five names are updated in `views.json` and appended to `state.json` `all_views`. The first run copies `views.json` to `views.json.bak`.

Thumbnails (`viewImage`) are not written here. Open each view in the MDV UI to generate them.

## How it works

1. **`infer_roles`** picks the observation datasource (usually `cells`), then matches column suffixes to roles: cell type, tissue, disease, treatment, sex, sample, QC, scores, UMAP pairs. RNA / protein links come from `rows_as_columns` subgroups.
2. **Markers** (optional) scan the sparse RNA matrix in `datafile.h5` for cluster-vs-rest top-20 genes and factor-varying genes. Results are cached as `cluster_markers.json` and a `cluster_markers` table datasource. Re-runs reuse the cache unless `--recompute-markers` is set.
3. **Feature UMAPs** on view 3 use unique top-2 genes per cluster (first cluster wins, cap 24). If markers were skipped, a canonical lineage panel is used instead.
4. Gene / protein charts use MDV wrappers `{subgroup}|{name}({subgroup})|{index}`.

## How to use it

After `mdvtools` is installed (repo checkout or `pip install -e python`):

```bash
mdvtools create-default-views /path/to/mdv_project
```

Skip the h5 marker scan (layout only):

```bash
mdvtools create-default-views /path/to/mdv_project --skip-markers
```

Recompute markers even when `cluster_markers.json` exists:

```bash
mdvtools create-default-views /path/to/mdv_project --recompute-markers
```

Module form (same flags as before):

```bash
python -m mdvtools.template_views --project /path/to/mdv_project
```

Gene panels, field-name patterns, and thresholds live in `defaults.py` and are imported by `create_default_views.py`.

Python API:

```python
from mdvtools.template_views import create_default_views, infer_roles

summary = create_default_views("/path/to/mdv_project", skip_markers=True)
```

## Tests

```bash
./venv/bin/python -m pytest python/mdvtools/tests/test_template_views.py
```
