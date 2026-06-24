# SpatialData.js integration notes (MDV)

Living log for MDV ↔ `@spatialdata/*` integration. Architecture and phased roadmap:

- [MDV integration](https://taylor-ccb-group.github.io/SpatialData.js/docs/vis/mdv-integration)
- [Headless viewer](https://taylor-ccb-group.github.io/SpatialData.js/docs/vis/headless-viewer)
- [MDV release checklist](https://taylor-ccb-group.github.io/SpatialData.js/docs/vis/mdv-release-checklist)

## State model

MDV owns `chart.config.renderStack` (MobX observable). The spatial chart composes rendering from `useSpatialCanvasRendererFromLayerInputs` + `SpatialViewer` (not the `SpatialCanvasViewer` black box) so MDV can publish the **Image Layer Registry** and wire **App Viv Extensions**. Tooltips remain MDV-owned (`renderTooltip={false}`); do not re-expose upstream internal tooltip paths that are not yet correct for MDV.

Chart root config owns chart/view concerns such as region, MDV overlay settings, and `viewState`. Render-stack entries own layer implementation concerns such as image channel settings (`entry.props.channels`). Viv channel/image settings are not part of the SpatialData chart root config.

The main integration concern is the MDV ↔ SpatialData.js render-stack boundary, not the local file count. MDV currently keeps a small adapter layer around that boundary to preserve object identity for spatial layer configs, cache cloned host deck layers, and force viewer refreshes without re-entering expensive geometry loads on cosmetic edits.

For this PR, keep that adapter local to MDV and make its invariants explicit before proposing SpatialData.js API changes. A follow-up SpatialData.js design opportunity is a render-stack adapter/hook that accepts a mutable stack plus a version token and returns identity-stable layer inputs and host deck layers.

Call the local boundary the **Render Stack Adapter**. It owns conversion from MDV's MobX-backed `config.renderStack` into SpatialData.js viewer inputs, including layer-config identity preservation, host deck layer clone caching, and explicit version-token refreshes. Viewer components should pass stack state into this adapter rather than coordinating cache invalidation, host fingerprints, and MobX observation directly.

Expose the adapter primarily as a React hook, backed by small pure helpers for testable transforms. The hook owns refs, cache lifetimes, memo dependencies, and host-layer fingerprinting; the chart component consumes viewer-ready `layers`, `layerOrder`, and `deckLayers`.

Keep the Render Stack Adapter read-side only. Stack editing remains in the render-stack control layer (`useRenderStackMutation`, `useRenderStackEntry`, insert/remove/reorder/default helpers), so config mutation behavior stays distinct from viewer input adaptation and performance caching.

Keep layer panels type-specific for this PR. A later pass can consolidate repeated panel controls, but this pass focuses on state adapters and helper boundaries around render-stack control and viewer input adaptation.

Organize `src/react/spatialdata` helpers around these responsibilities:

- **Render Stack Adapter** — read-side hook for SpatialData.js viewer inputs; owns layer-input cache, host fingerprinting, host clone caching, and generation-based refreshes.
- **Render Stack Control** — MobX-facing mutation hooks plus pure edit operations; owns patch/insert/remove/reorder behavior for `config.renderStack`.
- **Render Stack Defaults** — stack creation, normalization, available spatial entries, and default layer props. Prefer "defaults" over "seed" in names because it describes the domain role more directly.
- **Image Layer Runtime Bridge** — runtime-only helpers (`image_layer_runtime.ts`): `channelId`-keyed histogram stats, viewer UI array sync, tone → `vivLayerProps` mapping. Persistence is `useLayerChannelState` (upstream), not a local MobX bridge.

Leave orthogonal helpers separate: host overlay IDs, view-state bridging, spatial feature tooltip formatting, and table association placeholders.

Merge the current entry-state and mutation helpers into `render_stack_control.ts`, with pure stack operations first and MobX/React hooks second. This keeps all `config.renderStack` editing behavior under one named responsibility.

Keep `spatialEntryAsLayerConfig` with the adapter/read-side conversion code, not with mutation/control. It projects a stack entry into SpatialData.js `LayerConfig`; control hooks may use it for UI convenience, but it is not part of stack editing semantics.

Move host-layer fingerprinting and clone caching into the Render Stack Adapter. `host_overlay_ids.ts` remains the small vocabulary module for MDV overlay IDs, while the adapter owns the performance-sensitive mechanics that turn host stack entries into viewer-ready deck layers.

Move `createMdvHostLayerResolver` into the Render Stack Adapter too, so the chart component imports one module for the render-stack read path. Keep `host_overlay_ids.ts` purely declarative: MDV overlay IDs, labels, and host-layer ID conversion helpers.

Protect this refactor with a small number of focused Vitest unit tests for adapter/control/defaults behavior, especially object identity and cache invalidation. These tests should create useful friction against accidental performance regressions without exhaustively locking down implementation details. Prefer tests under `src/tests/react/spatialdata/`; defer Playwright coverage unless the user-visible dialog behavior changes.

Keep the adapter mostly MobX-neutral: the caller passes the mutable stack plus an explicit version token, and the adapter returns identity-stable viewer inputs. The token may remain named `renderStackGeneration`; the important point is that it exists to avoid garbage churn from replacing/cloning render-stack objects on cosmetic edits. Mutation/control remains the MobX-aware layer. The token should be clearly marked as an intentional, narrow refresh signal, not as a general state-management pattern.

Add a short code comment next to `renderStackGeneration` explaining that it is an intentional version token for in-place render-stack edits and exists to refresh adapter outputs without replacing stack objects or causing layer/data churn.

Do not leave compatibility re-export files for old helper names. These helpers are branch-local implementation details, so update imports to the new module boundaries and delete dead files rather than preserving aliases.

Remove the unused `refreshRenderStackShell` helper. Replacing only the stack shell can be revisited later if it proves cleaner without reintroducing layer/data churn, but this pass uses the explicit narrow refresh token.

Possible follow-up: consider publishing MDV to npm so SpatialData.js can exercise real MDV integration points without local worktree/link setup. This is an integration enabler, not part of the adapter refactor.

## Avivatorish comparison

MDV currently carries a local `src/react/components/avivatorish` implementation and also depends on `@spatialdata/avivatorish`. Treat the local copy as the integration shim for now, not as a desired long-term fork.

The two implementations are close enough to share vocabulary: `VivProvider`, `createVivStores`, `channelsStore`, `viewerStore`, `imageSettingsStore`, `useImage`, `useLoader`, `useMetadata`, and channel selection/state helpers. Both use Zustand-style stores under React providers, and both model channel state as parallel arrays (`colors`, `contrastLimits`, `domains`, `selections`, `channelsVisible`, etc.).

The important differences for MDV:

- MDV's local version is already wired into MDV chart/view concerns: MobX chart config, chart-link view state, OME-TIFF upload/viewer paths, and the existing Viv scatter/image chart controls.
- MDV's local channel stats path stores sampled `raster` data in channel state. The old channel histogram UI uses that raster data to draw histograms and brush contrast ranges.
- `@spatialdata/avivatorish` ships `useLayerChannelState`, `mergeLayerChannelState`, and `getSingleSelectionStats({ includeRaster: true })`.
- SpatialData.js exposes `useImageLayerContext`, `ImageLayerContextProvider`, and Viv extension passthrough on the composed renderer path.
- SpatialData.js image rendering is OME-Zarr oriented. MDV still has OME-TIFF paths on the legacy Viv chart, outside SpatialData image layers.
- MDV's local avivatorish copy remains the UI shim for histogram/tone controls on legacy `VivMdvReact`. The SpatialData chart uses upstream persistence hooks plus a thin runtime bridge.

Histogram brush controls are sensitive to state ownership. Keep the brush value controlled by one state source only: the channel store value that will be rendered (`contrastLimits` for image channels, or the equivalent legend range elsewhere). Avoid adding a local/debounced mirror that also writes through MobX or another persistence layer, because d3 brush movement, React re-render, and persistence rehydration can otherwise chase each other and cause drag jumps. If a bridge persists brush changes into another model, skip self-echo rehydration when the incoming persisted value matches the value just written.

Useful follow-up SpatialData.js changes (not blocking this PR):

- Array-of-structs `ChannelConfig` (deferred upstream ADR).
- Optional `onImageLayerRegistry` callback on `SpatialCanvasViewer` if MDV ever returns to the black-box viewer (MDV composes the renderer directly for now).

### Image layer panel pattern (MDV)

Dependencies: pin `@spatialdata/{core,layers,react,vis,avivatorish}` at **0.2.3**.

**Viewer (chart tree)**

1. `useRenderStackAdapter` → layer inputs (unchanged read-side cache).
2. `useSpatialCanvasRendererFromLayerInputs` + `SpatialViewer` — single load path; publish **Image Layer Registry** on `SpatialDataMdvReact` (`getImageLoadedDataByElementKey`, load state).
3. `vivImageExtensionResolver` / `vivImagePropsResolver` — `ColorPaletteExtension` + MDV `VivContrastExtension`; tone arrays from saved `entry.props.vivLayerProps`.
4. External tooltips only (`renderTooltip={false}`, MDV deck/feature tooltip portals). Do not wire upstream internal tooltip aggregation.

**Dialog (portal tree)**

[`ImageLayerPanel`](src/react/components/spatialLayers/ImageLayerPanel.tsx) per image stack entry:

1. `ImageLayerContextProvider` from chart registry + `elementKey` (provider colocated on panel, not dialog root).
2. `useImageLayerContext(elementKey)` — loader, defaults, `channelNames`, `selectionAxisSizes`. If `undefined`, show loading placeholder (no OME-metadata fallback).
3. `useLayerChannelState` — persistence to `entry.props.channels` via `patchLayer` / echo-safe hydrate.
4. Per-panel `VivProvider` — zustand projection for histogram/runtime only (`domains`, `raster`, loading flags).
5. **Spatial Image Panel Context** — exposes hook write API to `VivChannelList`; persisted edits (`colors`, `contrastLimits`, add/remove) go through the hook, not `channelsStore`.
6. [`image_layer_runtime.ts`](src/react/spatialdata/image_layer_runtime.ts) — `channelId`-keyed stats cache; rebuild index-aligned zustand arrays when hook state changes; tone sliders persist via `patchLayer({ vivLayerProps })`.

Delete `image_layer_channel_bridge.ts` and `use_image_layer_panel_defaults.ts` (replaced by upstream hook + registry).

**Legacy `VivMdvReact`** — unchanged; `ColorChannelComponents` branches on presence of Spatial Image Panel Context.

### Viv extension passthrough (landed upstream @0.2.3)

SpatialData.js `ChannelConfig` has an open TODO for channel-related extension props ([`packages/vis/src/SpatialCanvas/types.ts`](https://github.com/Taylor-CCB-Group/SpatialData.js/blob/main/packages/vis/src/SpatialCanvas/types.ts)). `VivSpatialViewer` must pass `extensions`, `brightness`, and `contrast` through `detailView.getLayers({ props })` without spreading `layer.props` afterward (drops non-enumerable extension props).

**Design split (same as host overlays):**

- **Serializable** on image stack entry: `channels` (core Viv fields via `useLayerChannelState`); `vivLayerProps` for tone (`brightness[]`, `contrast[]`) and other extension props.
- **Runtime** on composed renderer: `vivImageExtensionResolver` / `vivImagePropsResolver` — MDV supplies `LayerExtension` instances.

**Other upstream (done or deferred):** `useLayerChannelState`, `useImageLayerContext`, `includeRaster` stats — **0.2.3**. Array-of-structs `ChannelConfig` deferred.

MDV root `config.viv` is **not** an upstream concern.

## Initial PR scope

| Commit stage | Status |
|--------------|--------|
| 1 | `@spatialdata/*@0.2.3` deps + CONTEXT + this doc |
| 2 | `SpatialDataMdvRegionReact` — composed renderer + host overlays + Image Layer Registry |
| 3 | Layer dialog — `renderStack` list + dnd-kit reorder |
| 4 | Visibility + opacity per stack entry |
| 5+ | Image layer panel — registry, `useLayerChannelState`, runtime bridge, `vivLayerProps` tone |

Implemented under `src/react/spatialdata/` and `src/react/components/SpatialData*`.

## Deferred (follow-up PR)

- Table-driven shape colouring (`fillColorByColumn`, `spatial_table_association`)
- `@spatialdata/avivatorish` zarr loader delegation (MDV keeps OME-TIFF local)
- Playwright fixture test

## Divergence from prototype branch

`codex/spatialdata.js_first_pass` used deprecated `layers`/`layerOrder` and `SpatialLayerStackConfig`. This worktree uses library `RenderStack` only.
