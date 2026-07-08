# 05 — Async / GPU / WASM filtering and untangling `useFilteredIndices`

> Making Dimension filtering async (WebGPU/WASM-backed), cleaning up the `useFilteredIndices`
> hook, and building filter-graph primitives — with reusable GPU code from `tgpu-htj2k`.

## The one-line diagnosis

**Dimension filter *evaluation* is the main synchronous main-thread hotspot.** The *worker-backed
aggregations* derived from a filter (histograms, category counts, box/violin, density contours,
sort order, and the compact passing-index list) are **already off-thread in workers reading
SharedArrayBuffers**. So the biggest win from "make filtering async" = move the
predicate-over-all-rows step to a worker/GPU/WASM.

Two important caveats so this isn't overstated: (i) the **chart-scope predicate** and **cross-chart
ownership** work in the `useFilteredIndices` cluster still run **synchronously on the main thread**
(see §"what crept in" items 3 and 5) — they are derived-from-filter work that is *not* off-thread;
(ii) the Dimension ref-count bookkeeping is main-thread too. The worker/SAB architecture is a good
substrate for GPU/WASM, but "everything derived is already off-thread" applies to the Dimension
aggregation workers, not to the React-side scoping/ownership layer.

## Current filter model (grounded)

- A `DataStore` owns one global filter byte-array (`filterArray` over a `SharedArrayBuffer`,
  [DataStore.js:74](../../../src/datastore/DataStore.js)); a row passes iff `filterArray[row] === 0`.
  The global array is actually a **per-row reference count of excluding dimensions**; `filterSize`
  tracks passing rows.
- Each active filter is a `Dimension` ([Dimension.js](../../../src/datastore/Dimension.js)) owning
  its own local byte array with states `0/1/2/3` (0 none, 1 local, 2 background, 3 both).
  `_applyStateTransition` (Dimension.js:18) does incremental ref-count bookkeeping — only when a
  row's "is excluded" boolean flips does it touch the parent count.
- **Filter evaluation is 100% synchronous, main-thread, full O(size) scans**: `filterPredicate`,
  `filterRange/Square/Poly`, `filterCategories`, `filterValueset`, `filterCatCol`,
  `setBackgroundFilter`. Dimension.js:167 even `console.warn`s when a filter exceeds 100ms.
- Comment at Dimension.js:45 already muses about "a new strategy for filter evaluation that doesn't
  attempt to use this class" and passing in `filteredIndices` to avoid full scans.

### Sync vs async today

| Operation | Today |
|---|---|
| Predicate/range/poly/category **filter evaluation** | **SYNC main thread**, full scan |
| Background filter apply/clear, `removeAllFilters`, `reFilterOnDataChanged` | **SYNC main thread** |
| Compact passing-index list (`getFilteredIndices` → `filteredIndexWorker.ts`) | ASYNC (worker), promise-cached |
| Filtered/unfiltered histograms, category counts, dotplot, box/violin, density, sort, heatmap | ASYNC (workers over SABs) |

All workers **read** filter state from SABs and only aggregate; none writes the filter. Filter
state is **byte-per-row in shared memory** — a flat, transferable, structured-clone-free buffer,
ideal for GPU/WASM.

## `useFilteredIndices` — what it is and what crept in

`useFilteredIndices()` ([src/react/hooks.ts:477](../../../src/react/hooks.ts)) returns the
`Uint32Array` of currently-passing row indices. Precisely: the **shared, worker-computed compact
index list** is produced by `useSimplerFilteredIndices` (hooks.ts:493) — that is the piece that
turns the 35M-row byte array into a cached ~1M-element list shared across charts. `useFilteredIndices`
wraps it and **layers the synchronous chart-scope predicate on top** (hooks.ts:386), so the value it
returns is not the raw shared list but a per-chart-scoped subset. The performance win (deck.gl /
Splatter / DeckScatter only draw passing points) comes from the `useSimplerFilteredIndices` compact
list and is genuinely load-bearing for the spatial rendering path; the extra scoping layer is part
of what should be separated out (below).

But the single conceptual job ("give me passing indices") has accreted ~9 distinct
responsibilities across `hooks.ts`, `filterOwnership.ts`, `categoryFilterUtils.ts`, and a divergent
copy in `scatter_state.ts`:

1. DataStore→React event bridging (`useFilteredIndicesRefreshVersion`, hooks.ts:361).
2. Async worker orchestration + observable promise-cache mirroring (`useSimplerFilteredIndices`,
   hooks.ts:493).
3. Per-chart "chart-scope" category/background filtering — an **entirely separate**, in-memory,
   main-thread filter layered on top of the Dimension system, with its own async column loading +
   readiness state machine (`useChartScopeFilterPredicate`, hooks.ts:386).
4. Column-load readiness gating (`loadedFilterKey` vs `filterKey`) — the source of surprising
   zero-length arrays.
5. Cross-chart filter-ownership math (`useOwnedFilteredIndices` + `filterOwnership.ts:59`) —
   full-array scans on the main thread.
6. Single-category sub-selection deliberately bypassing Dimensions for speed
   (`useCategoryFilterIndices`, hooks.ts:566 — comment explains a Dimension "always passes through
   all rows, which this doesn't").
7. Raw byte-array passthrough with a synthetic dependency (`useFilterArray`).
8. Background-filter state-byte semantics leaking into React helpers (`localFilter[row] === 2`).
9. A parallel, divergent `scatter_state.ts` implementation with different config assumptions.

## `tgpu-htj2k` primitives worth reusing

Sibling repo [`tgpu-htj2k`](https://github.com/xinaesthete/tgpu-htj2k) (public; may be renamed) —
WebGPU (TypeGPU/TGSL) + Rust→WASM. No barrel export; import by path. Paths below are relative to
that repo. Directly reusable:

- **`getDevice()`** (`src/gpu/device.ts`) — cached GPU device; in-browser backed by `navigator.gpu`.
- **The layout-bound compute-pipeline idiom** (`src/gpu/spatial/nnDistance.ts`):
  `tgpu.bindGroupLayout` with `params`(uniform) / `pts`(readonly storage) / `outb`(mutable
  storage), `tgpu.computeFn` authored in TS (`"use gpu"` TGSL), workgroup size 64, pooled/grown
  storage buffers. **This is the template for a per-row GPU filter kernel** — one thread per row,
  write 0/1 to an output buffer.
- **`splatDensityGpu(xs, ys, {weights})`** (`src/gpu/spatial/splatDensity.ts`) — rasterizes a
  **weighted point cloud into an r32float density grid** via additive blending. This *is* a GPU 2D
  histogram / KDE. Set `weights[i] = passing ? 1 : 0` and you get a **live filtered density field
  with zero extra marshalling** — the x/y column buffers and the filter byte array are already in
  shared memory. Natural backend for the Splatter chart and a "filter heatmap."
- Rust→WASM pipeline (`build:wasm --target web`) demonstrated by `htj2kCore.ts` — the packaging
  path a WASM filter kernel could reuse (not a filtering primitive itself).

## Opportunities

### A. Async Dimension filter evaluation

1. **Row-predicate GPU kernel.** Column data is already flat typed arrays over SAB, and
   `Dimension.filterArray` is a flat `Uint8Array` over SAB. A compute kernel modeled on
   `nnDistance.ts` takes `(columnBuffer, params) → filterArray` writing 0/1 per row. Insertion
   point: a new evaluation strategy behind `Dimension.filter(method, …)`. Range/square/valueset are
   embarrassingly parallel; polygon (`filterPoly`) ray-casting maps to a per-row loop over verts.
2. **Keep ref-count semantics off the GPU.** The tricky part is `parent.filterArray` ref-count +
   `filterSize`. Options: (a) GPU produces the new *local* byte array; a parallel reduction
   computes the delta and new `filterSize`; the main-thread merge becomes a bulk diff; or (b)
   **recompute the global array as a parallel sum of dimension arrays** (states 1/3 = excluded) —
   this simplifies away the incremental ref-counting entirely.
3. **WASM fallback.** For no-WebGPU environments, a Rust→WASM predicate pass over the same SAB in a
   worker gives a big constant-factor win off the main thread.
4. **Reify filter kinds as serializable data.** `filterPredicate`'s JS-closure predicate can't
   cross into WASM/GPU — restructure the `filter(method, columns, args)` dispatch so `args` are
   serializable (range bounds, category-id set, polygon verts) rather than opaque functions. This
   is a prerequisite for both GPU and WASM paths.
5. **Reuse existing plumbing.** `getFilteredIndices` / `filteredIndexWorker` already prove the
   SAB-in/SAB-out worker round-trip + MobX promise-cache pattern (the single shared
   `_filteredIndexWorker` per store, DataStore.js:854); a GPU-backed `filterWorker` can mirror it.
   Both filter workers carry a `//todo atomics` note — lock-free double-buffering is the goal.

### B. Filter-histogram / filter-density graphs

1. **GPU density straight from filter state.** `splatDensityGpu(xs, ys, {weights: passingMask})`
   yields a live filtered density field with no extra data marshalling — backend for the Splatter
   chart and filter-heatmap overlays.
2. **GPU/parallel histograms.** `rawHistogramWorker.ts` already has a TODO to accept
   `filteredIndices`; a GPU histogram (atomic-add into bins, or the render-to-accumulate trick from
   `splatDensity`) over `col.buffer` gated by the filter mask gives instant filtered histograms —
   replacing the flaky `getBinsAsync` (RangeDimension.js:156, self-described "wasted hours … gave
   up").
3. **A unified "filter graph" reduction.** Today each readout (bins, cats, box, density) has a
   bespoke worker + a bespoke `getX(callback,…)` on a Dimension subclass. A single GPU/WASM
   reduction over `(columnBuffer, filterMask, reductionSpec)` could collapse
   `binWorker`/`catWorker`/`catColWorker`/`boxPlotWorker` into one parameterized primitive —
   shrinking the subsystem and making filtered graphs cheap.

### C. Untangling `useFilteredIndices`

1. **Separate the three fused layers:** (i) store→React signal
   (`useFilteredIndicesRefreshVersion`, keep); (ii) global passing indices
   (`useSimplerFilteredIndices`, keep as the async/worker source of truth); (iii) per-chart scoping
   (`useChartScopeFilterPredicate` + ownership math) — this is really a **local filter** concept
   and should be modeled as a proper (possibly GPU-backed) Dimension or a first-class "chart-local
   filter," not an in-memory `.filter(predicate)` bolt-on. The top-of-file comments (hooks.ts:471)
   already argue for this generalization.
2. **Kill the divergent `scatter_state.ts` copy** by making the general local-filter mechanism the
   single implementation.
3. **Move whole-array scans out of React** (`filterOwnership.ts:59`, `useOwnedFilteredIndices`)
   into the same worker/GPU path that produces the compact list — return "my visible rows" and
   "externally filtered" as buffers rather than recomputing per render on the main thread.
4. **Reify chart-scope filters as Dimensions** so filter graphs (B) automatically reflect them —
   today `background_filter`/`category_filters` live only in React config and are invisible to the
   worker aggregations, which is why only dotplot/scatter honor background filters.

## Relationship to the rest

This theme is **largely orthogonal** to the spatial-tables work, but it shares the same
performance concern (the 35M-cell / ~1M-visible rendering path) and the same architectural grain
(flat SABs, worker offload). The cleanup (C) can precede the GPU work (A/B). Reifying filter kinds
as serializable data (A4) and chart-scope filters as Dimensions (C4) are the two changes that
unlock the most downstream value and should come first.

## Suggested order within this theme

1. **Reify filter args as serializable data** (A4) + **reify chart-scope filters as Dimensions**
   (C4) — prerequisites that also clean up the model.
2. **Untangle `useFilteredIndices`** into the three layers (C1–C3); kill the `scatter_state.ts`
   duplicate.
3. **GPU filter-density graphs** via `splatDensityGpu` (B1) — highest visible payoff, self-contained.
4. **GPU/WASM predicate evaluation** (A1–A3) with the global-array-as-parallel-sum simplification.
5. **Unified filter-graph reduction** (B3) — the big consolidation, once the primitives are proven.
