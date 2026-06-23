# SpatialData.js integration notes (MDV)

Living log for MDV ↔ `@spatialdata/*` integration. Architecture and phased roadmap:

- [MDV integration](https://taylor-ccb-group.github.io/SpatialData.js/docs/vis/mdv-integration)
- [Headless viewer](https://taylor-ccb-group.github.io/SpatialData.js/docs/vis/headless-viewer)
- [MDV release checklist](https://taylor-ccb-group.github.io/SpatialData.js/docs/vis/mdv-release-checklist)

## State model

MDV owns `chart.config.renderStack` (MobX observable). `SpatialCanvasViewer` is a controlled renderer: plain `renderStack` + `hostLayerResolver` at the React boundary. No parallel stack maps or periodic whole-stack snapshots.

Chart root config owns chart/view concerns such as region, MDV overlay settings, and `viewState`. Render-stack entries own layer implementation concerns such as image channel settings (`entry.props.channels`). Viv channel/image settings are not part of the SpatialData chart root config.

The main integration concern is the MDV ↔ SpatialData.js render-stack boundary, not the local file count. MDV currently keeps a small adapter layer around that boundary to preserve object identity for spatial layer configs, cache cloned host deck layers, and force viewer refreshes without re-entering expensive geometry loads on cosmetic edits.

For this PR, keep that adapter local to MDV and make its invariants explicit before proposing SpatialData.js API changes. A follow-up SpatialData.js design opportunity is a render-stack adapter/hook that accepts a mutable stack plus a version token and returns identity-stable layer inputs and host deck layers.

Call the local boundary the **Render Stack Adapter**. It owns conversion from MDV's MobX-backed `config.renderStack` into SpatialData.js viewer inputs, including layer-config identity preservation, host deck layer clone caching, and explicit version-token refreshes. Viewer components should pass stack state into this adapter rather than coordinating cache invalidation, host fingerprints, and MobX observation directly.

Expose the adapter primarily as a React hook, backed by small pure helpers for testable transforms. The hook owns refs, cache lifetimes, memo dependencies, and host-layer fingerprinting; the chart component consumes viewer-ready `layers`, `layerOrder`, and `deckLayers`.

Keep the Render Stack Adapter read-side only. Stack editing remains in the render-stack control layer (`useRenderStackMutation`, `useRenderStackEntry`, insert/remove/reorder/default helpers), so config mutation behavior stays distinct from viewer input adaptation and performance caching.

Keep layer panels type-specific for this PR. A later pass can consolidate repeated panel controls, but this pass focuses on state adapters and helper boundaries around render-stack control and viewer input adaptation.

Organize `src/react/spatialdata` helpers around three responsibilities:

- **Render Stack Adapter** — read-side hook for SpatialData.js viewer inputs; owns layer-input cache, host fingerprinting, host clone caching, and generation-based refreshes.
- **Render Stack Control** — MobX-facing mutation hooks plus pure edit operations; owns patch/insert/remove/reorder behavior for `config.renderStack`.
- **Render Stack Defaults** — stack creation, normalization, available spatial entries, and default layer props. Prefer "defaults" over "seed" in names because it describes the domain role more directly.

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
- `@spatialdata/avivatorish` exposes channel stats helpers, but the current declared stats contract returns `domains` and `contrastLimits`, not the sampled raster data needed by MDV's histogram component.
- SpatialData.js `SpatialCanvas` image loading exposes loaded image defaults (`colors`, `contrastLimits`, `channelsVisible`, `selections`) for layer panels, but not a public histogram/raster-stat surface for channel-control UI.
- SpatialData.js image rendering is OME-Zarr oriented. Its image renderer notes that SpatialData image support uses `loadOmeZarr`; MDV still has OME-TIFF paths that are outside that renderer's current scope.
- MDV's local version has awkward Zustand/MobX mixing because Viv state is runtime UI state while `chart.config.renderStack` is serialized MobX state. This is tolerable as a bridge, but it should not leak Viv root config back into SpatialData chart config.

For the current SpatialData chart, image channel controls may reuse the old MDV channel components only for state that can be faithfully backed by `renderStack.entries[].props.channels`. Until SpatialData.js exposes real channel histogram/raster samples, MDV may render the existing histogram brush over empty raster data so users can still edit `contrastLimits`; replace that placeholder with a SpatialData.js-provided stats API when one exists.

Histogram brush controls are sensitive to state ownership. Keep the brush value controlled by one state source only: the channel store value that will be rendered (`contrastLimits` for image channels, or the equivalent legend range elsewhere). Avoid adding a local/debounced mirror that also writes through MobX or another persistence layer, because d3 brush movement, React re-render, and persistence rehydration can otherwise chase each other and cause drag jumps. If a bridge persists brush changes into another model, skip self-echo rehydration when the incoming persisted value matches the value just written.

Useful SpatialData.js changes before MDV can consider replacing the local `avivatorish` copy:

- Expose a public channel-control adapter or hook that takes a layer-local image channel config and returns editable channel state plus a persistence callback, without requiring root Viv chart config.
- Expose image channel histogram/stat data, or a lazy per-selection stats API, alongside loaded image defaults.
- Make the supported image source boundary explicit: either support OME-TIFF in the shared package, or keep OME-TIFF intentionally MDV-local and separate from SpatialData image layers.
- Clarify whether the shared package supports Viv/deck.gl extensions needed by MDV (`VivContrastExtension`, color palette/colormap extensions, 2D/3D behavior), and expose extension injection where host apps need it.
- Keep the channel ownership model layer-local: image channel settings belong to image layer configs/render-stack entries, not to the SpatialData chart root.

## Initial PR scope

| Commit stage | Status |
|--------------|--------|
| 1 | `@spatialdata/*@next` deps + CONTEXT + this doc |
| 2 | `SpatialDataMdvRegionReact` chart — `SpatialCanvasViewer` + host overlays |
| 3 | Layer dialog — `renderStack` list + dnd-kit reorder |
| 4 | Visibility + opacity per stack entry |
| 5+ | Type panels; shapes/labels/points added via dialog only |

Implemented under `src/react/spatialdata/` and `src/react/components/SpatialData*`.

## Deferred (follow-up PR)

- Table-driven shape colouring (`fillColorByColumn`, `spatial_table_association`)
- `@spatialdata/avivatorish` zarr loader delegation (MDV keeps OME-TIFF local)
- Playwright fixture test

## Divergence from prototype branch

`codex/spatialdata.js_first_pass` used deprecated `layers`/`layerOrder` and `SpatialLayerStackConfig`. This worktree uses library `RenderStack` only.
