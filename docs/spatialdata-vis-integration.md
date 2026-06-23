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
