# 00 — First-class table ↔ spatial-element association (the core deliverable)

> **This is the headline user-facing feature.** The other docs (01–05) are enablers or
> adjacent work. The tangible win: a spatial **element** (shapes / labels) that is annotated by a
> **table** gets a proper, first-class link to that table's rows — so the geometry is **colored,
> filtered, highlighted, and tooltipped by table columns**, and filtering flows both ways between
> the geometry and every other chart on that table.

## What "first-class" means (the proven reference: the points path)

MDV already has a first-class table↔geometry association — for the **points** representation. The
SpatialData region scatter renders a table's rows as deck.gl points and is first-class because
**three concerns share one per-row index space, all rooted on the `DataStore`**
([src/react/scatter_state.ts](../../../src/react/scatter_state.ts)):

- **Shared filter mask.** `data` fed to the layer is `ownerVisibleRows` — a `Uint32Array` of
  surviving row indices from `useOwnedFilteredIndices()`
  ([hooks.ts:520](../../../src/react/hooks.ts)). Only passing rows draw; brushing *any* chart
  updates `ds.filterArray` and the geometry redraws. Lasso on the spatial view feeds back via
  `spatial_context.tsx` → `rangeDimension.filterPoly` → `ds.filterArray`.
- **Shared color pipeline.** `getFillColor: colorBy`, where
  `colorBy = chart.getColorFunction(col, true)` returns `(rowIndex) => [r,g,b]`
  ([BaseChart.ts:683](../../../src/charts/BaseChart.ts) → [DataStore.js:1524](../../../src/datastore/DataStore.js)).
  `updateTriggers.getFillColor: colorBy` recolors by function-reference identity. A legend is
  produced as a side effect.
- **Shared highlight.** Deck picking `onClick → dataHighlighted(data[index])`
  ([selectionHooks.ts:133](../../../src/react/selectionHooks.ts)); highlight rings render from
  `useHighlightedIndices()`.

**The goal is to give shapes/labels geometry the exact same three properties**, keyed on the same
`DataStore` row-index space — so a shape and a point representing the same cell filter, colour, and
select together.

## The two halves that must meet

### MDV side — per **row** (already exists)
`ds.filterArray` (shared filter) · `chart.getColorFunction(col, true)` (row → colour) ·
`dataHighlighted(rowIndex)` (shared highlight) · tooltip via `config.tooltip.column`. All indexed
by **DataStore row**.

### Library side — per **feature** (already exists, in SpatialData.ts)
The association layer lives in
`SpatialData.ts/packages/core/src/tableAssociations.ts` and is element-agnostic:

- `loadAssociatedTableFeatureRows({spatialData, kind, key, extraColumnNames})` resolves the
  annotating table via `getAssociatedTable` (matching the element key against the table's `region`
  array), then builds `rowIndexByFeatureId: Map<instanceKeyString, absoluteObsRow>`, filtered by
  the `region_key` column. Critically `TableElement.loadObsIndex()` returns the **`instance_key`**
  column as row ids (not `obs.index`) — matching Python spatialdata's join-by-instance semantics.
- `loadFeatureRowIndexByFeatureIndex({..., featureIds})` → **`Int32Array`** answers "rendered
  feature *i* → table row" (`-1` = unmatched).
- For shapes, `ShapesElement.loadRenderData()` returns `ShapesRenderData` with
  `featureIds: string[]` and **`rowIndexByFeatureIndex: Int32Array` pre-filled**. Geometry loads
  from `shapes/<key>/shapes.parquet` (parquet-wasm → arrow, WKB-decoded) in parquet row order.
- The library's shapes deck layers (`PolygonLayer` / `ScatterplotLayer`, **not** `GeoJsonLayer`)
  already carry `{featureId, featureIndex, rowIndex}` on every datum, expose a per-feature
  colouring surface `featureState.fillColorByFeatureId` / `hiddenFeatureIds` / `fadedFeatureIds`,
  and emit pick events that include `rowIndex`
  (`SpatialData.ts/packages/layers/src/shapesLayer.ts`).

### The bridge
`rowIndexByFeatureIndex[i]` maps rendered feature *i* → DataStore row. Invert it (via `featureIds`)
to go row → featureId when handing colours to the library's `featureState` (which is
**featureId-keyed, not row-keyed**).

### Alignment guarantees (from the library)
- **Instance-key join** is the primary contract: feature id (shape GeoDataFrame index / label
  pixel value) ↔ the table row whose `instance_key` equals it, restricted by `region_key`.
- **Positional fallback (shapes only)** fires *only* when feature ids are exactly `"0".."n-1"` (a
  synthesized parquet RangeIndex) and row counts match — then feature *i* ↔ the *i*-th filtered
  table row. Fragile; assumes parquet order == filtered-table order.
- Unmatched features get `-1` and are surfaced explicitly, not guessed.

## Why nothing connects today

The two halves exist but are not wired:

- `table_association.ts` ([src/react/spatialdata/table_association.ts](../../../src/react/spatialdata/table_association.ts))
  is a **type stub** (`none | ambiguous | resolved`) with **no resolver** — `NO_TABLE_ASSOCIATION`
  is hardcoded at both panel call sites.
- The `fillColorByColumn` and `tooltipFields` UI in
  [ShapesLayerPanel.tsx:138](../../../src/react/components/spatialLayers/ShapesLayerPanel.tsx) is
  **commented out** ("not working yet"); the `LabelsLayerPanel` is **disabled**.
- Shapes render with a **static fill**, ignore `ds.filterArray`, and are not picked into the
  DataStore. Their tooltip is library-driven, not MDV column-driven.
- MDV supplies **no per-feature accessors** to the library's shapes layer today — only static
  `fillColor` / `strokeColor` props via the render-stack `LayerConfig`.

## The Python contract gap (why this pulls toward the JS-read path)

The association joins on **`instance_key`** (or positional row order). The current `convert-spatial`
path **destroys both**:

- `convert-spatial` reads `instance_key` via `get_table_keys` but **discards it** — `_instance_key`
  at [conversion.py:373](../../../python/mdvtools/spatial/conversion.py) and `:639`. The only
  per-row id that survives is `mdv_cell_id`, derived from the **merge-suffixed** obs index
  ([conversions.py:242](../../../python/mdvtools/conversions.py)).
- The single-`anndata.concat` merge **row-stacks all tables** ([02-convert-spatial-tables.md](02-convert-spatial-tables.md)),
  scrambling row order — so even positional alignment can't be trusted.

**Consequence:** the clean association resolves against the **untouched zarr store**, where
SpatialData.ts's association layer already computes the feature→row map. So the headline feature is
most naturally delivered on the **JS-read path** — a spatial chart backed by the SpatialData object
directly ([01-js-dataloader.md](01-js-dataloader.md)) plus the shared store cache
([03-shared-spatial-contexts.md](03-shared-spatial-contexts.md)). The alternative is to fix the
converter (theme 2) to **preserve `instance_key` and `region_key` as columns and stop scrambling
row order** (one-datasource-per-table / merge-by-region) so the h5 path can join too.

> This is the single most important strategic conclusion of the session: **the core deliverable
> and the JS DataLoader / converter rework are the same project, not separate ones.** You cannot
> have first-class association on top of the current lossy h5 flattening.

## Rendering options

Both need the **prerequisite**: an association resolver (fill in `table_association.ts`) that maps
an element (`kind`, `elementKey`, coordinate system) → the MDV `DataStore` for its table, and
builds the feature→row `Int32Array` (via `ShapesElement.loadRenderData()` /
`loadFeatureRowIndexByFeatureIndex`).

### Option A — drive the library's `featureState` (reuse library geometry *and* rendering)
Compute, from MDV per-row data, the featureId-keyed inputs the library layer already consumes:
- `fillColorByFeatureId[featureId] = colorBy(rowIndexByFeatureIndex[i])` using
  `chart.getColorFunction(col, true)` — **the same colour pipeline as points**.
- `hiddenFeatureIds` / `fadedFeatureIds` from `ds.filterArray` (a feature is hidden/greyed iff its
  row is filtered out) — **the same filter mask as points**.
- wire the library's shape pick event (`rowIndex`) → `dataHighlighted` and `config.tooltip.column`.

*Pro:* reuses the library's geometry loading *and* deck layers; least new rendering code.
*Con / must-verify:* requires `@spatialdata/vis`'s renderer (`useSpatialCanvasRendererFromLayerInputs`)
to accept per-layer `featureState` (or equivalent colour/hidden inputs) through the render-stack
`LayerConfig` / renderer inputs. The `@spatialdata/*@0.2.5` packages are **not installed in this
worktree** — confirm this injection seam against installed `node_modules` before committing to A.

### Option B — MDV host deck layer using library-loaded geometry (recommended for parity)
Call `ShapesElement.loadRenderData()` to get geometry + `rowIndexByFeatureIndex` (so **no WKB /
parquet reimplementation**), then build an **MDV-owned** `PolygonLayer` / `ScatterplotLayer` as a
**host overlay** — MDV already injects host deck layers (`deck:scatter`, gates, selection) into the
spatial canvas via `createMdvHostLayerResolver`
([render_stack_adapter.ts](../../../src/react/spatialdata/render_stack_adapter.ts),
[host_overlay_ids.ts](../../../src/react/spatialdata/host_overlay_ids.ts)). The overlay mirrors the
points path verbatim:
- `data` = features whose mapped row is in `ownerVisibleRows` (shared filter),
- `getFillColor: (f) => colorBy(rowIndexByFeatureIndex[f.featureIndex])` with
  `updateTriggers.getFillColor: colorBy` (shared colour),
- `pickable`, `onClick → dataHighlighted(row)` + highlight stroke from `useHighlightedIndices()`
  (shared highlight),
- tooltip via `config.tooltip.column` over the mapped row.

*Pro:* genuinely first-class — one `filterArray`, one `getColorFunction`, one `dataHighlighted`
shared with points; no dependence on the library's (unverified) colour-injection surface.
*Con:* a new host-layer category (generalize the fixed `deck:*` id enum); MDV owns a second shapes
renderer alongside the library's.

**Recommendation:** prototype **A** first (cheapest if the seam exists); fall back to **B** for
guaranteed parity. Either way the colour machinery to reuse is `chart.getColorFunction(col, true)`
→ `(row) => [r,g,b]`, wired by function reference — the exact seam the points layer uses.

### Labels are harder (raster-only)
Labels have **no vector geometry** — only a segmentation raster; "features" are integer pixel
values (label ids == instance_key), resolved to rows only at pick time
(`tooltipRowIndexByFeatureId.get(String(labelId))`, instance-key join, no positional fallback). So:
- **Colour by table column** = build an `instanceKey → [r,g,b]` LUT from `getColorFunction` and have
  the label layer recolour the raster by that map (check whether `@spatialdata/vis`'s labels layer
  accepts a colour map; `useLayerData` has a labels path).
- **Filter** = a hidden-label mask (paint filtered rows to background/transparent).
- You **cannot enumerate all label features** without scanning the raster.
Do labels **after** shapes; treat colour-by-LUT as the MVP and filtering as a follow-up.

## Bidirectional filtering

- **Table → geometry:** `ds.filterArray` change → recompute hidden/visible features → redraw
  (Option A `hiddenFeatureIds`, or Option B `data = ownerVisibleRows`). Free once the bridge exists.
- **Geometry → table:** lasso/box on the spatial view → point-in-polygon over shape
  centroids (or the picked feature set) → `rangeDimension.filterPoly` → `ds.filterArray`, exactly
  as the points path already does via `spatial_context.tsx`.

## Phasing

0. **Association resolver** — fill `table_association.ts`: map `elementKey` → the table's DataSource
   (via region metadata), build the feature→row `Int32Array`. Everything else depends on this.
1. **Colour shapes by a table column** — reuse `getColorFunction`; drive Option A `featureState` or
   Option B host layer. Re-enable the `fillColorByColumn` UI in `ShapesLayerPanel`.
2. **Filter shapes by the table** — `filterArray` → hidden/faded features; redraw.
3. **Pick + highlight + tooltip shapes** → `dataHighlighted`, `config.tooltip.column`.
4. **Geometry → table lasso selection** (bidirectional filtering).
5. **Labels**: colour-by-LUT, then hidden-label filtering.

## Open questions / risks

- **Which DataSource is "the table" for an element?** The resolver must map `elementKey` → region →
  the MDV datasource name. Clean on the JS path (the store's `getAssociatedTable`); needs
  converter-preserved region metadata on the h5 path.
- **The converter identity gap** (above) — decide JS-read vs converter-fix as the delivery route.
- **`@spatialdata/vis` `featureState` injection seam** — verify against installed packages
  (worktree lacks them) before betting on Option A.
- **Multi-region tables**: the library's `getAssociatedTable` uses only the *first* match; a
  broad "matched zero → use all rows" fallback can over-associate for shared tables.
- **Positional-alignment fragility** — safe only for `0..n-1` parquet indices with matching row
  counts; another reason to prefer real `instance_key`s end to end.
- **Perf**: `loadFeatureRowIndexByFeatureIndex` re-scans the whole obs table per call; cache it in
  the project-scoped store cache (theme 3) for interactive use.
- **Labels raster limits** — no feature enumeration; colour/filter via LUT + masking only.
