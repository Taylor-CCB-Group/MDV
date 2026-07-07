# 03 — Project-scoped shared contexts (image / store cache)

> Why every spatial chart currently re-opens the zarr store and re-fetches images, and the safe
> path to a project-scoped cache. This is the enabling substrate for theme 1's DataLoader, and
> (per your framing) "one of the safer things to address."

## The goal

Two spatial charts (or a chart and its layer dialog) pointing at the same SpatialData path should
share **one opened store** and **one image cache**, instead of each independently calling
`readZarr` and re-loading pixels.

## Why it's duplicated today

**Every MDV chart and dialog is an independent React island.** `BaseReactChart.mountReact()` →
`createMdvPortal()` does a fresh `createRoot` per chart/dialog
([src/react/components/BaseReactChart.tsx:76](../../../src/react/components/BaseReactChart.tsx),
[src/react/react_utils.tsx:98](../../../src/react/react_utils.tsx)). There is no shared root; the
class docstring flags "a single root with portals" as unfinished future work
([BaseReactChart.tsx:34](../../../src/react/components/BaseReactChart.tsx)).

**But a project-scoped provider pattern already works.** `createMdvPortal` wraps every island in:

```
<QueryClientProvider client={queryClient}>        // module singleton, react_utils.tsx:78
  <ChartManagerProvider chartManager={window.mdv.chartManager}>
    <ProjectProvider …>
      {component}
```

The Provider *elements* are duplicated per island, but their *values* are **singletons**, so the
same context value crosses island boundaries. This is the template.

### Three cache layers, all per-island today

| Layer | What | Where | Shared? |
|---|---|---|---|
| (a) `SpatialData` store (zarr) | `useMemo(() => readZarr(source), [source])` | upstream `SpatialDataProvider`, mounted at [SpatialDataMDVReactComponent.tsx:297](../../../src/react/components/SpatialDataMDVReactComponent.tsx) **and** [SpatialLayerDialogReactWrapper.tsx:30](../../../src/react/components/SpatialLayerDialogReactWrapper.tsx) | **No** — memo is per provider instance; chart tree + dialog tree each open their own |
| (b) store-instance caches | `parquetTableCache`, `obsIndices`, `varIndices`, … | live *on* the `SpatialData` object (`@spatialdata/core`) | No — each `readZarr` yields a distinct instance |
| (c) Viv image loader / loaded pixels | OME-Zarr multiscale data | owned by `useSpatialCanvasRendererFromLayerInputs` ([…Component.tsx:124](../../../src/react/components/SpatialDataMDVReactComponent.tsx)) | No — one hook instance per chart viewer |

The **Image Layer Registry**
([src/react/spatialdata/image_layer_registry.ts](../../../src/react/spatialdata/image_layer_registry.ts))
bridges (c) from a chart's viewer tree to *its own* layer dialog — it is a **within-one-chart**
bridge, not cross-chart. Chart B's dialog cannot see chart A's registry.

Confirmed by grep: there is **no** `spatialStore`, `SpatialStoreProvider`, `spatialDataCache`, or
shared `readZarr` wrapper in `src/`. Two charts on the same `region.spatial.file` re-open and
re-fetch.

## The provider/context inventory (what NOT to touch)

| Context | Holds | Scope | Verdict |
|---|---|---|---|
| `QueryClientProvider` / `ChartManagerProvider` / `ProjectProvider` | singletons | app/project | **precedent to follow** |
| `VivProvider` / `vivStores` (zustand) | runtime channel/UI state (histograms, brush) | per-chart (and a fresh copy per dialog panel) | **do not share** — intentionally per-chart; for the SpatialData chart it's not even the source of truth (canonical state is MobX `renderStack`) |
| `SpatialDataProvider` (`spatialDataPromise`) | `readZarr(source)` | per provider instance | **the one to dedupe** |
| `SpatialAnnotationProvider` | MDV scatter/gate/selection deck layers | per-chart | keep per-chart (chart-specific overlays) |
| `ImageLayerContextProvider` / `SpatialImagePanelContext` | registry callbacks / per-panel write API | per image panel | keep panel-local |

## Recommended path

### Step 1 — project-scoped SpatialData store cache (safe, high payoff)

- New module `src/react/spatialdata/spatial_store_cache.ts`: a `Map<string, Promise<SpatialData>>`
  keyed by resolved store URL (+ `selection` if ever used). One instance per project — a module
  singleton mirroring `queryClient`, or hung off `chartManager` / `ProjectContext`.
- Add a thin `SpatialStoreProvider` inside `createMdvPortal`
  ([react_utils.tsx:108](../../../src/react/react_utils.tsx)), alongside `ProjectProvider`,
  exposing the cache.
- Replace the two direct `<SpatialDataProvider source>` usages with a wrapper that resolves the
  store from the shared cache (or wrap upstream `SpatialDataProvider` so its `readZarr` is deduped
  by URL). Result: N charts + dialogs on one path → **1 `readZarr`, 1 `SpatialData` instance,
  shared `parquetTableCache` etc.** (layers (a) and (b) shared) with minimal surface area.

This alone requires **no single-root refactor** — the shared value crosses islands via the
singleton, exactly like `queryClient`.

### Step 2 — share the Viv image loader / loaded pixels (larger)

Today loaded-image data is owned by one renderer hook and exposed via the per-chart registry. To
reuse pixel/multiscale loads across charts you need a project-scoped image-loader cache keyed by
`(store, elementKey, selection)`, then have each chart's renderer consult it. This likely needs an
upstream `@spatialdata/vis` affordance (the integration doc already flags a "registry accessor on
composed renderer" gap, [docs/spatialdata-vis-integration.md:159](../../spatialdata-vis-integration.md)).
Higher risk — do after step 1.

### Step 3 — connect to the DataLoader (theme 1)

The project-scoped `SpatialData` store from step 1 is the **same object** a JS DataLoader reads
table columns from. Placing the cache on ChartManager/Project scope makes it reachable from both
React (via context) and non-React DataLoader code (via `window.mdv.chartManager`), the way
`ProjectContext` is already consumed both ways.

## Risks

- **React island boundaries** — share via a **singleton value**, never via React tree ancestry.
- **MobX vs zustand** — canonical layer state is MobX (`config.renderStack`,
  `renderStackGeneration`), runtime channel UI is zustand. Keep the store cache a **plain
  singleton/Map**, *not* MobX-observable, or you risk re-render churn in the identity-sensitive
  render-stack adapter ([render_stack_adapter.ts](../../../src/react/spatialdata/render_stack_adapter.ts)).
- **Lifetime / eviction** — a project-scoped cache outlives individual charts; define eviction on
  project unload / ChartManager teardown, or accept session lifetime like `queryClient`. Watch for
  leaks of large zarr/image data.
- **StrictMode double-invoke** — `createMdvPortal` renders under `<StrictMode>`; populate the cache
  in the memo/getter (idempotent by URL key), not in an effect.
- **Registry teardown timing** — the Image Layer Registry is set/cleared in an effect tied to
  `renderStackGeneration`. If image loading moves to a shared owner (step 2), revisit who
  publishes/tears down the registry so the dialog panel still resolves.

## Why this is the safe first move

The provider stack already demonstrates working project-scoped singletons crossing island
boundaries. A store cache is **additive**, follows that exact precedent, is orthogonal to the
MobX/zustand machinery, and (step 1) touches only two provider mount sites plus one new module —
no single-React-root refactor, no change to the perf-critical render-stack adapter. And it unblocks
theme 1.
