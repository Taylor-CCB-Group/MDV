# 01 — A JS (spatialdata.js / anndata.js) DataLoader with h5 write routing

> Read [README.md](README.md) first for the four load-bearing facts. This doc details the
> data-layer contract, what a zarr read loader must produce, and how it reconciles with the
> existing h5 write path.

## The goal

Add a `DataLoader` variant whose `.function` reads column data directly in the browser from a
SpatialData / AnnData **zarr** store (via `spatialdata.js` / `anndata.js`), instead of POSTing
to the Flask `/get_data` endpoint. Some columns remain **h5-backed and writable** — especially
anything the user edits — so the loader must know which columns to read from zarr and which to
defer to the server, and writes must continue to land in `datafile.h5`.

## Current data layer (grounded)

### The DataLoader contract

`getDataLoader(isStaticFolder, datasources, views, url)`
([src/dataloaders/DataLoaderUtil.ts:85](../../../src/dataloaders/DataLoaderUtil.ts)) returns:

```ts
type DataLoader = {
  function:        (columns: ColumnInfo[], dataSource: string, size: number)
                     => Promise<Array<{ field: string; data: SharedArrayBuffer }>>,  // real contract
  viewLoader:      (view: string) => Promise<any>,
  rowDataLoader:   (dataSource: string, index: string) => Promise<any>,
  binaryDataLoader:(dataSource: string, name: string) => Promise<ArrayBufferLike>,
}
```

> The declared type of `.function` (`Promise<ArrayBufferLike>`) is wrong/loose; the runtime
> contract is a **per-column list of `{field, data}`**. Worth fixing while here.

Today there are two variants, selected by `isStatic`
([projectRuntime.ts:61](../../../src/modules/projectRuntime.ts)):
- **Server** (`isStatic=false`): `getArrayBufferDataLoader('{root}/get_data')` — one POST per
  batch, returns a single concatenated buffer split by `processArrayBuffer`
  ([DataLoaders.ts:197](../../../src/dataloaders/DataLoaders.ts)).
- **Static** (`isStatic=true`): `getLocalCompressedBinaryDataLoader` — per-column HTTP Range
  requests against `{ds}.gz` with a `{ds}.json` byte-offset index
  ([DataLoaders.ts:238](../../../src/dataloaders/DataLoaders.ts)).

### The invariants a loader must honor

`processArrayBuffer` ([DataLoaders.ts:106](../../../src/dataloaders/DataLoaders.ts)) is the
canonical byte→typed-array spec. Any loader must produce, per column, a `SharedArrayBuffer`
laid out as:

| datatype | TypedArray | notes |
|---|---|---|
| `integer` / `int32` | `Int32Array` | length `size` |
| `double` | **`Float32Array`** | MDV stores double as **float32** — downcast float64 |
| `text` | `Uint8Array` | length `size`; values are **indexes into `column.values`** |
| `text16` | `Uint16Array` | length `size`; indexes into `values` |
| `unique` | `Uint8Array` | fixed width `stringLength`; length `size*stringLength` |
| `multitext` | `Uint16Array` | width `stringLength`, `65535` sentinel for empty |

Plus:
- **Column-major**: one contiguous buffer per column, exactly `size` rows. `size` is fixed per
  datasource from `datasources.json` ([DataStore.js:74](../../../src/datastore/DataStore.js)).
- **SharedArrayBuffer**, not `ArrayBuffer` — data is posted to workers without copying. (Requires
  cross-origin isolation / COOP+COEP; verify headers if reading from a CDN.)
- **Categorical columns store integer indexes into a `values: string[]` array**, and `values`
  must be present in metadata and **stable across read and write** or colors/filters break.
- **Sparse matrix columns** expand to a dense `Float32Array(size)` pre-filled with `NaN`.
- The **4-byte alignment ordering** (numeric columns first) only matters if you build a single
  concatenated buffer. **Return one SAB per column and it is a non-issue** — recommended.

### Where loaded bytes land

`ChartManager._loadColumnData` ([ChartManager.js:1495](../../../src/charts/ChartManager.js))
calls `this.dataLoader(columns, dataSource, size)`, then for each returned `{field, data}` calls
`dataStore.setColumnData(field, data)` ([DataStore.js:1262](../../../src/datastore/DataStore.js)),
which accepts a `SharedArrayBuffer` directly. `columnsWithData.push(field)` marks it loaded.
**None of this cares where the bytes came from.**

### The write path (unchanged by any read loader)

Edits: `getState()` → `_getUpdatedColumns` (materializes typed arrays back to JS arrays,
[ChartManager.js:1101](../../../src/charts/ChartManager.js)) → `state_saved` listener →
`POST /save_state` → `save_state` → `set_column_with_raw_data`
([mdvproject.py:849](../../../python/mdvtools/mdvproject.py)) → h5 `create_dataset`. The write
target is **always `datafile.h5`**. Writability is whole-project (`MDVProject.writable`,
[mdvproject.py:139](../../../python/mdvtools/mdvproject.py)); there is no server-side per-column
`editable` flag (the client sends one but the server ignores it).

## What the JS loader must implement

A new factory, e.g. `getZarrDataLoader(datasources, storeDescriptor)`, returning the `DataLoader`
shape. Its `.function`, given `(columns, dataSource, size)`:

1. Resolve the datasource's zarr store (ideally from the **project-scoped store cache** — see
   [03-shared-spatial-contexts.md](03-shared-spatial-contexts.md)).
2. For each requested column, read `size` values keyed by `column.field` (or
   `subgroup`/`sgindex` for matrix-backed rows-as-columns).
3. Encode to the exact typed-array layout above, wrap in `SharedArrayBuffer`, return
   `{field, data}[]`.
4. Provide `viewLoader` / `rowDataLoader` / `binaryDataLoader` — these can **delegate to the
   existing server loaders**; only column data comes from zarr.

### What the JS libraries give us (and don't)

From the read of `anndata.js` and `SpatialData.ts` (both **read-only**, both on `zarrita`):

- **Whole-column reads are cheap and the natural mode.** `SpatialData.ts`
  `TableElement.loadObsColumns([name])` → `VAnnDataSource._loadColumn`
  (`SpatialData.ts/packages/core/src/models/VAnnDataSource.ts:86`) returns plain JS
  `TableColumnData` with **categoricals already decoded to `string[]`** and per-column promise
  caching. This is the **most directly reusable** path — closest analogue to MDV `get_data`.
- `anndata.js` (`obs.get(name)` + `get(handle, [null])`) returns zarrita `Chunk`s; better for
  **X-matrix** slicing (dense and sparse). Sparse densifies a whole **major-axis** slice — pick
  CSR vs CSC to match the hot access direction (cell-wise vs gene-wise).
- **No arbitrary row-index gather** in either library. If `get_data` is asked for a filtered row
  subset, load the whole column and index in JS, or drop to zarrita for chunk-aware gather.
- Column discovery: `TableElement.getObsColumnNames()` (sync, from the parsed tree). `anndata.js`
  has no "list columns" method.
- **Version skew risk:** `anndata.js` pins `zarrita@0.5.1`; `SpatialData.ts` targets `0.7.x`.
  Reconcile before sharing a bundle. `anndata.js` is `0.0.x` — expect API churn.
- MDV already depends on `@spatialdata/*` (`^0.8.0` as of the table-association work) and
  `zarrextra` — so the packages are on hand. `zarrextra` is no longer image-only: the association
  reads obs columns through the store, though **table column data still loads over h5**, which is
  the gap this theme is about.

### The one net-new client concept: `backing`

The client has **no model today for "read from zarr, write to h5."** Introduce column metadata,
e.g. `backing: "h5" | "zarr"` (or `readSource`/`writeTarget`), threaded through:
- `datasources.json` schema (`DataSourceSchema`, validated in
  [DataStore.js:29](../../../src/datastore/DataStore.js)),
- `DataColumn` type ([charts.d.ts:75](../../../src/charts/charts.d.ts)),
- `getColumnInfo` ([DataStore.js:1876](../../../src/datastore/DataStore.js)),
- `_getUpdatedColumns.getMd` ([ChartManager.js:1155](../../../src/charts/ChartManager.js)) so the
  save payload can carry it if the server ever needs it.

The loader consults `backing` to **split the column list**: zarr for `backing:"zarr"`, HTTP
`/get_data` for the rest; merge the `{field,data}[]` results. `_loadColumnData` passes the whole
list and just consumes the merged result, so the fan-out is internal to the loader.

## Recommended design: a hybrid read loader

- `.function` receives the full column list, partitions by `backing`, reads zarr columns via
  `SpatialData.ts`'s `loadObsColumns` (whole-column), reads the rest via the existing `/get_data`
  POST, encodes both to SABs, and merges.
- **Staleness rule:** a user edit writes the column to h5 and upserts its metadata into
  `datasources.json`. The simplest rule: **once a column exists as an h5 dataset, prefer h5.**
  `/get_configs` already returns live datasource metadata after a save, so an h5-shadowed column
  can be flagged there and the loader routes it to `/get_data` on the next load.
- **Server changes are optional.** `/save_state`→`set_column_with_raw_data` already writes any
  column to h5. You only touch the server if you want edits mirrored back into zarr, or want the
  server to track the h5-shadow set for read routing.

### Encoding pitfalls to plan for

- **`double` → float32** downcast: watch precision on columns that round-trip zarr → edit → h5.
- **Categorical `values` ordering** must be identical on read and write; `anndata.js` /
  `SpatialData.ts` eagerly decode codes→strings, but MDV wants **codes + a `values` dictionary**
  — consider reading the `codes`/`categories` zarr arrays directly to preserve the compact form
  rather than re-deriving `values` from decoded strings.
- **`unique` fixed width**: zarr strings are variable-length; compute `stringLength = max` and pad
  to match `_convertColumn` / `processArrayBuffer`.
- **Sparse X**: MDV hands `DataStore` a dense `NaN`-filled `size`-length SAB even for sparse —
  memory-heavy for large gene matrices; a zarr loader can stream but must still densify per column.

## Phasing

1. **Read-only proof of concept.** New `getZarrDataLoader`; branch in `getDataLoader` keyed on a
   per-datasource `store` descriptor (stub type `ExperimentalZarrStore` exists at
   [charts.d.ts:161](../../../src/charts/charts.d.ts)). All columns from zarr; no writes yet.
   Validate numeric + categorical columns against a known project.
2. **Hybrid routing.** Add `backing` metadata; fan out zarr vs `/get_data`. Confirm an edited
   column round-trips to h5 and reloads from h5.
3. **Store-cache integration.** Read from the project-scoped `SpatialData` cache (theme 3) so the
   loader and the spatial charts share one opened store.
4. **Matrix / rows-as-columns.** Serve gene-expression subgroups from zarr X (choose CSR/CSC).

## Open questions / risks

- **Source of truth after an edit** (zarr vs h5) — the "prefer h5 if dataset exists" rule needs a
  reliable signal; decide whether the server or the client owns the shadow set.
- **SharedArrayBuffer cross-origin isolation** when zarr is served from a different origin.
- **Row-index gather** absence — acceptable if `get_data` stays whole-column (it is today), but
  revisit if async filtering (theme 5) wants to fetch only passing rows.
- **Version alignment** of zarrita across `anndata.js` / `SpatialData.ts` / MDV.
- **Greenfield**: `@spatialdata/*` is wired only to image rendering today — no prior art for the
  column-data path in MDV.
