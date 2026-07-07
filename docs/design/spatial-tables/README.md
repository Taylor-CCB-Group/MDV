# SpatialData tables & runtime-data roadmap — research and strategy

> **Status:** exploratory research, not decisions. This folder is a design-space map to
> support scoping conversations. Nothing here is committed; the ADRs in `docs/adr/` are the
> place for decisions that are hard to reverse.
>
> Written from a deep read of the codebase (frontend `src/`, backend `python/mdvtools/`) plus
> the sibling libraries `SpatialData.ts`, `anndata.js`, and `codecs/tgpu-htj2k`. Every claim is
> anchored to `file:line` so you can verify quickly.

## The core deliverable

**A spatial element (shapes / labels) annotated by a table gets a proper first-class association
to that table's rows — so the geometry is coloured, filtered, highlighted, and tooltipped by
table columns, and filtering flows both ways between the geometry and every other chart on that
table.** This is the tangible user-facing win. → **[00-table-element-association.md](00-table-element-association.md)**

MDV already does this for the **points** representation of a table (the region scatter shares one
`DataStore` row-index space for filter + colour + highlight). The gap is that **shapes/labels
geometry does not** — the association is a stub (`TableAssociation` type with no resolver), the
`fillColorByColumn` UI is commented out, and shapes render with a static fill, ignore the filter,
and are not picked into the DataStore.

The other themes below are **enablers or adjacent work** for that deliverable — most importantly,
the clean association joins on `instance_key` against the **untouched zarr store**, which the
current h5-conversion path destroys. So the core feature and the read-path / converter work are
one project (see [00](00-table-element-association.md) → "the Python contract gap").

## The enabling / adjacent themes

1. **A new `DataLoader`** that reads tabular data via `spatialdata.js` / `anndata.js`
   (`get_data` from zarr), while still routing writes to the h5 datafile. *(Enabler: the
   association resolves natively on the JS-read store.)*
   → [01-js-dataloader.md](01-js-dataloader.md)
2. **Flexible table handling in `convert-spatial`** — stop force-merging incompatible tables
   (grid + cells) into one datasource, and **preserve `instance_key` / `region_key`** so the h5
   path can also join. *(Enabler / alternative delivery route.)*
   → [02-convert-spatial-tables.md](02-convert-spatial-tables.md)
3. **Project-scoped shared contexts** so the image/zarr cache works across spatial charts
   instead of every chart re-opening the store. *(Enabler: one opened store for the association
   + image loads.)* → [03-shared-spatial-contexts.md](03-shared-spatial-contexts.md)
4. **Rethinking views / DataSources / ChartManager** — a datasource picker in "Add Chart",
   cross-datasource shared layout, a new view kind, and swapping the ChartManager at runtime.
   *(Adjacent: how a table datasource and its element share screen space.)*
   → [04-views-datasources-chartmanager.md](04-views-datasources-chartmanager.md)
5. **Async / GPU / WASM filtering** and untangling `useFilteredIndices`, with filter-graph
   primitives from `tgpu-htj2k`. *(Adjacent: the same 35M-cell filter path the geometry shares.)*
   → [05-dimension-async-filtering.md](05-dimension-async-filtering.md)

## The four load-bearing facts (read these first)

These recur across every theme and largely determine what is cheap vs. expensive.

- **The DataLoader seam is single and clean.** `getDataLoader(isStatic, datasources, views,
  url)` ([src/dataloaders/DataLoaderUtil.ts:85](../../../src/dataloaders/DataLoaderUtil.ts))
  returns one object; the load-bearing member is `.function(columns, dataSource, size) →
  Promise<{field, data: SharedArrayBuffer}[]>`. Everything downstream
  (`ChartManager._loadColumnData`, `DataStore.setColumnData`) is source-agnostic. A new loader
  plugs in here and nothing else needs to know.

- **Writes always go to h5, regardless of read source.** Edits flow
  `getState()` → `updatedColumns` → `POST /save_state` → `set_column_with_raw_data` → h5
  `create_dataset` ([python/mdvtools/mdvproject.py:849](../../../python/mdvtools/mdvproject.py)).
  Writability is a **whole-project** filesystem check, not per-column. So a zarr *read* loader
  coexists with the h5 *write* path for free — the only new concept needed is a `backing`
  hint on column metadata for **read routing** (zarr vs `/get_data`), not for permissions.

- **Almost everything is keyed by DataSource *name*.** `views.json` is `viewName → {dataSources,
  initialCharts}` keyed by ds name; `addChart(dsName, config)` takes a name; a chart's ds
  binding is *only* which `initialCharts[dsName]` array it sits in — the chart config never
  names its datasource ([src/charts/ViewManager.ts:9](../../../src/charts/ViewManager.ts),
  [src/charts/ChartManager.js:2247](../../../src/charts/ChartManager.js)). This is why the
  datasource-picker is trivial, the runtime ChartManager swap is plausible, and "datasources
  derived from spatial objects" can bind to existing views — as long as names/columns match.

- **Project-scoped singletons already cross React islands.** Every chart/dialog is its own
  `createRoot`, but `createMdvPortal` wraps each island in `QueryClientProvider` /
  `ChartManagerProvider` / `ProjectProvider` whose *values* are singletons
  ([src/react/react_utils.tsx:98](../../../src/react/react_utils.tsx)). A shared spatial store
  cache follows this exact precedent — no single-root refactor needed.

## How the themes connect

```
        ╔═════════════════════════════════════════════════════════════╗
        ║  CORE DELIVERABLE (00): shapes/labels ↔ table association    ║
        ║  geometry coloured / filtered / highlighted BY table columns ║
        ║  reuse the points path: one DataStore row-index space        ║
        ╚════════════▲══════════════════▲══════════════════▲══════════╝
                     │ joins on instance_key against the untouched store
        ┌────────────┴───────┐  ┌───────┴────────────┐  ┌──┴───────────────┐
        │ JS DataLoader (1)  │  │ convert-spatial (2)│  │ shared store /    │
        │ zarr read; assoc.  │  │ preserve instance_ │  │ image cache (3)   │
        │ works natively     │  │ key; per-table dss │  │ one opened store  │
        └──────┬─────────────┘  └────────────────────┘  └───────────────────┘
               │  (delivery route: JS-read OR converter-fix)
               ▼
        ┌─────────────────────────────────────────────────────────────┐
        │  Views / DataSources / ChartManager (4)                      │
        │  datasource picker · cross-DS layout · new view kind ·       │
        │  runtime ChartManager swap (rebind same views to new data)   │
        └─────────────────────────────────────────────────────────────┘

        ┌─────────────────────────────────────────────────────────────┐
        │  Dimension filtering (5): async/GPU/WASM + useFilteredIndices│
        │  the SAME filter mask the geometry shares (35M-cell path)    │
        └─────────────────────────────────────────────────────────────┘
```

Observations worth calling out:

- **The association is the goal; 1/2/3 are how you make it clean.** The clean join needs real
  `instance_key`s against the original store. That is native on the **JS-read path** (1) with the
  **shared store cache** (3), or achievable on the **h5 path** only if the **converter** (2)
  preserves `instance_key`/`region_key` and stops scrambling row order.
- **Theme 3's store cache is also theme 1's data source** and the association layer's cache — the
  project-scoped `SpatialData` store that dedupes image loads is the same object the DataLoader
  reads columns from and `loadFeatureRowIndexByFeatureIndex` scans. Build it first.
- **The "radical" runtime-derived-datasources vision is themes 1+2+4 composed** — read tables
  directly (1) or emit clean per-table datasources (2), synthesize DataSource configs, and swap
  the ChartManager while keeping views (4). Each is independently useful.

## Risk / reward at a glance

| Change | Blast radius | Reward | Confidence |
|---|---|---|---|
| Association resolver + colour shapes by column (00) | New resolver + shapes render seam | **Core deliverable** | Medium — library provides the join; needs the right store |
| Filter/highlight/pick shapes via DataStore (00) | Host layer or library `featureState` | **Core deliverable** | Medium — mirrors the points path |
| Labels colour/filter (00) | Raster LUT + masking | High | Low-Medium — raster-only, no feature enumeration |
| Datasource picker in Add Chart (4) | Dialog only | Medium (UX) | High — `DatasourceDropdown` already exists, unused |
| Project-scoped SpatialData store cache (3) | 1 new module + 2 mount sites | High (perf, enables 1) | High — follows `queryClient` precedent |
| convert-spatial table grouping (2) | Python converter, offline | High (correctness) | High — localized to one merge block |
| JS zarr read DataLoader (1) | 1 new loader + `backing` metadata | High (no re-convert) | Medium — encoding/format edge cases |
| Async/GPU filter evaluation (5) | Dimension subsystem | High (35M-cell path) | Medium — ref-count semantics are fiddly |
| `useFilteredIndices` untangling (5) | React hooks + consumers | Medium (maintainability) | Medium — ~9 entangled concerns |
| Cross-DS shared layout (4) | Panes, gridstack, getState, view schema | Medium (UX) | Medium — structural but clean |
| New view kind bypassing ChartManager (4) | `_init`/`getState`/`hasUnsavedChanges` | High (flexibility) | Low-Medium — many implicit assumptions |
| Runtime ChartManager swap (4) | Singleton + portals teardown | High (runtime vision) | Low-Medium — singleton coupling |

## Suggested sequencing

A dependency-respecting order that front-loads the substrate the **core deliverable** needs, then
builds the association feature on top:

1. **Project-scoped SpatialData store cache** (theme 3, step 1). Safe, high payoff; the substrate
   the DataLoader *and* the association layer both read from. Independent of everything else.
2. **Decide the delivery route for identity** — JS-read (theme 1) vs converter-fix (theme 2). The
   association joins on `instance_key`; pick where clean `instance_key`s come from. A minimal
   read-only JS table loader (1) against the shared store (1←3) is the shortest path to a clean
   join without re-converting data; the converter fix (2) is the durable h5-path answer.
3. **Association resolver + colour shapes by a table column** (theme 00, phases 0–1). Fill
   `table_association.ts`; reuse `getColorFunction`. The first visible slice of the deliverable.
4. **Filter + highlight + pick shapes via the DataStore** (theme 00, phases 2–3). Full parity with
   the points path — the heart of "first-class."
5. **Bidirectional lasso selection + labels** (theme 00, phases 4–5). Labels last (raster-only).
6. **Adjacent, sequence by appetite:** datasource picker (4B), async filtering + `useFilteredIndices`
   cleanup (5), then the invasive view/ChartManager changes (4B/A/C) if the runtime-derived
   vision is confirmed.

## Per-theme documents

| # | Doc | Core question |
|---|---|---|
| **00** | **[00-table-element-association.md](00-table-element-association.md)** | **How do shapes/labels become coloured/filtered/highlighted by their table? (the deliverable)** |
| 1 | [01-js-dataloader.md](01-js-dataloader.md) | How does a zarr-backed read loader coexist with h5 writes? |
| 2 | [02-convert-spatial-tables.md](02-convert-spatial-tables.md) | How to group tables into datasources instead of one concat? |
| 3 | [03-shared-spatial-contexts.md](03-shared-spatial-contexts.md) | How to share the image/store cache across charts? |
| 4 | [04-views-datasources-chartmanager.md](04-views-datasources-chartmanager.md) | How to break per-DS segregation and swap data at runtime? |
| 5 | [05-dimension-async-filtering.md](05-dimension-async-filtering.md) | How to make filtering async and untangle `useFilteredIndices`? |
