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

## What is wired today

> Rewritten after the fact. The original text of this section — "why nothing connects today" —
> described a stub resolver, commented-out UI, and static fills. None of that is still true, so
> it has been replaced rather than annotated. Everything above this line was written before the
> work and still reads true.

`table_association.ts` has a real resolver, both panels are live, and geometry is coloured,
filtered and tooltipped by table columns. What did **not** get built is the return direction:
picking a shape does not reach `dataHighlighted`, and lasso on the spatial view does not select
shapes into the DataStore.

Colour turned out to have **two routes**, decided per column by one question — is the column in
the annotating table's `obs`? — which neither Option A nor Option B below anticipated:

- **In obs** → hand the column to the viewer as `fillColorByColumn`, with MDV's own palette
  attached (`fillColorSchemeFromDataStore`). The viewer reads the column, encodes it, and does its
  own load-window retention. This needed upstream work: `fillColorByColumn` could originally only
  cycle an index-ordered list, which cannot express "this category is this colour".
- **Not in obs** (gene scores, `mdv_cell_id`, anything computed at runtime) → MDV is the only
  party holding the values, so it colours each feature itself via `featureState.fillColorByFeatureId`
  — Option A as described.

Filtering is always `featureState`: cross-filtering is MDV's, and no column in obs can express it.
The full routing table is in [docs/spatialdata-vis-integration.md](../../spatialdata-vis-integration.md).

So Option A's seam was real, and the answer to "does `@spatialdata/vis` accept per-layer
`featureState`?" is yes — but the more interesting finding is that for a column the viewer *can*
read, handing over the column beats handing over a `Record<string, [r,g,b,a]>` with one entry per
cell on every filter change.

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
`LayerConfig` / renderer inputs.

> **Resolved: A, and the seam is real.** Verified against `@spatialdata/*@0.6.0` — the floor for
> this work, because the palette and domain fields it needs do not exist below it. See
> [What is wired today](#what-is-wired-today) for the obs/not-obs split that came out of it.

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

| | Phase | State |
|---|---|---|
| 0 | **Association resolver** — map `elementKey` → the table's DataSource (via region metadata), build the feature→row `Int32Array` | **Done** |
| 1 | **Colour shapes by a table column** — reuse `getColorFunction`; drive Option A `featureState` or Option B host layer | **Done** (Option A, plus the obs route) |
| 2 | **Filter shapes by the table** — `filterArray` → hidden/faded features; redraw | **Done** |
| 3 | **Pick + highlight + tooltip shapes** → `dataHighlighted`, `config.tooltip.column` | **Tooltip done; highlight open** — the library emits pick events carrying `rowIndex`, but nothing routes them into `dataHighlighted` |
| 4 | **Geometry → table lasso selection** (bidirectional filtering) | **Open** |
| 5 | **Labels**: colour-by-LUT, then hidden-label filtering | **Done** — labels take the same two colour routes as shapes |

Labels arriving early rather than last is worth noting: they came almost free once colour was
expressed as *a column the viewer resolves* rather than *a map MDV computes*, because the viewer
already knew how to read a column against a labels element. The doc's assumption that labels
would need a bespoke LUT path held only for the not-in-obs route.

## Columns that come from elsewhere (links, vars, computed)

**Any numeric column in MDV can accept a `RowsAsColsLink`** — that is the rule, and it is why the
column picker offers the "active link" tab wherever the parameter accepts numeric. A spatial layer's
fill colour is no exception; it now takes one (`field_spec_projection.ts`), so a layer can be
coloured by "whichever gene is selected over there" and follow that selection live.

Worth naming for future work: **the way `RowsAsColsLink` is used today is really "choose a `var`
from a table"** — the linked datasource is a gene/feature table, and the link picks a column of the
expression matrix. Reading it that way suggests two directions, neither of which changes the
current implementation:

- **Computed columns.** A column produced by an expression graph is the same shape of problem — a
  column identity that resolves late and can change while the chart is open. `mdvFieldSpecs` (the
  spec kept on MDV's side of the layer props, projected onto the viewer's concrete field names) is
  already the seam for that: a different kind of spec, same projection.
- **Cheaper `vars`.** Resolving a var column currently materialises it across the whole datasource.
  An annotated element often covers **far fewer rows than the datasource has** — one region of a
  multi-region table — so the values actually needed are a small slice. Scoping the fetch to the
  element's rows is the obvious win, and it wants the feature→row `Int32Array` this theme already
  builds. Nothing today is structured to take advantage of that.

## Open questions / risks

- ~~**Which DataSource is "the table" for an element?**~~ Answered on the JS path:
  `resolveAssociatedElementTable`, via the store's `getAssociatedTable`. The h5 path still needs
  converter-preserved region metadata.
- ~~**The converter identity gap**~~ — settled by delivery: JS-read. The converter gap (theme 2) is
  sidestepped, not closed, and still bites anything that wants the association on the h5 path.
- ~~**`@spatialdata/vis` `featureState` injection seam**~~ — verified against `0.6.0`; the seam
  exists and Option A was taken.
- **Multi-region tables**: still open, and now visible in the routing. `obsColumnNamesForElement`
  only answers when there is exactly **one** associated table, and returns `undefined` otherwise —
  which quietly sends every column of an ambiguously-associated element down the per-feature route
  rather than mis-colouring it. That is a safe default, not a solution.
- **Positional-alignment fragility** — safe only for `0..n-1` parquet indices with matching row
  counts; another reason to prefer real `instance_key`s end to end.
- **Perf**: `loadFeatureRowIndexByFeatureIndex` re-scans the whole obs table per call; cache it in
  the project-scoped store cache (theme 3) for interactive use.
- **Labels raster limits** — no feature enumeration; colour/filter via LUT + masking only.
