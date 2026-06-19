# SpatialData.js integration notes (MDV)

Living log for MDV ↔ `@spatialdata/*` integration. Architecture and phased roadmap:

- [MDV integration](https://taylor-ccb-group.github.io/SpatialData.js/docs/vis/mdv-integration)
- [Headless viewer](https://taylor-ccb-group.github.io/SpatialData.js/docs/vis/headless-viewer)
- [MDV release checklist](https://taylor-ccb-group.github.io/SpatialData.js/docs/vis/mdv-release-checklist)

## State model

MDV owns `chart.config.renderStack` (MobX observable). `SpatialCanvasViewer` is a controlled renderer: plain `renderStack` + `hostLayerResolver` at the React boundary. No parallel stack maps or periodic whole-stack snapshots.

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
