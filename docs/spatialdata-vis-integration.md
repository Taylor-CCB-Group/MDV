# SpatialData.js integration notes (MDV)

This document captures integration friction between MDV and `@spatialdata/*`, the
dual-tier design goals for SpatialData.js, and MDV's near-term workarounds.

## Design goals

SpatialData.js should serve two audiences:

| Tier | Audience | Needs |
|------|----------|-------|
| **Low-level** | Apps on `@spatialdata/core` + `@spatialdata/layers` with custom deck.gl composition | Stable layer types, `useSpatialCanvasRenderer`, explicit `layerOrder`, documented host overlay ids |
| **High-level** | Apps wanting turnkey `SpatialCanvas` | Exported headless panel primitives, store serialization, bundled renderer context |

MDV sits between both: it composes scatter, gates, and selection overlays on top of
spatial layers, and provides a richer layer-management dialog locally.

## Avivatorish triplication

Viv/image channel code exists in three places:

| Copy | Location | Loader routes |
|------|----------|---------------|
| MDV local | `src/react/components/avivatorish/` | OME-TIFF, multi-TIFF, OME-Zarr, JPEG2000 |
| Package | `@spatialdata/avivatorish` | OME-Zarr / SpatialData zarr only |
| Via vis | `@spatialdata/vis` → `@spatialdata/avivatorish` | Same as package |

**MDV shim strategy (near-term):**

1. Keep MDV local `avivatorish/` for `VivMDVReact` and OME-TIFF paths.
2. Use `LayerConfig.channels` for SpatialData zarr images in the layer dialog.
3. Share presentational channel UI via source-aware `ImageLayerPanel` adapters.
4. Optionally delegate zarr loader creation to `@spatialdata/avivatorish` from MDV
   `createLoader` in a follow-up.

**Longer-term:** restore OME-TIFF as optional loader plugins in `@spatialdata/avivatorish`;
contribute hardened MDV panel primitives upstream as headless building blocks.

## Layer stack (MDV)

`SpatialDataMdvReact` persists `config.spatialLayerStack`:

- Spatial layers (`image`, `shapes`, `points`, `labels`) keyed by layer config id
- Deck overlays keyed as `deck:grey_scatter`, `deck:scatter`, etc.
- Unified `stackOrder` drives `layerOrder` passed to `SpatialViewer`

Legacy GeoJSON `json` overlay is retired; region geometry uses spatial `shapes` layers.

## Table ↔ datastore association

Conversion merges all SpatialData tables into one MDV obs datasource with columns:

- `spatialdata_path` — zarr store basename
- `spatial_region` — `{store}_{element}`
- `table_name` — SpatialData table key
- `instance_key` — feature join key (from `uns.spatialdata_attrs`)

**Multi-table limitations:**

- `spatial_region` does not encode `table_name`
- Region metadata does not record table↔element bindings
- `background_filter` scopes by region only

**MDV happy-path inference** (`src/react/spatial_table_association.ts`):

1. Resolve store from `region.spatial.file`
2. Resolve element from shapes `elementKey` or parse `config.region`
3. Collect `table_name` candidates from filtered MDV rows
4. Collect candidates from SpatialData table attrs annotating the element
5. Single-table stores resolve immediately; intersection of candidates resolves;
   otherwise `ambiguous` or `none`

**Planned future:** per-table datasource or spatialdata-native dataloader; emit
`associated_table` in region metadata at conversion time.

## Upstream improvement backlog

1. Export headless panel primitives (`@spatialdata/vis/panels`)
2. Pluggable `createLoader` with optional OME-TIFF / multi-TIFF plugins
3. Shared `ChannelConfig` type + `useImageChannelState` hook
4. `ExternalLayerDescriptor` for host deck overlays in `layerOrder`
5. Replace naive `composeSpatialDeckLayers` concat with `composeLayers({ entries, layerOrder })`
6. Stable `orderId` on external layers; dev warnings for unknown `layerOrder` ids
7. `SpatialCanvasRendererProvider` bundling store + renderer hook outputs
8. `serializeStore()` / `hydrateStore()` on spatial canvas store
9. Composable viewer primitives instead of black-box `SpatialCanvasViewer`
10. `AssociatedTableInfo` on shapes/labels layer config

## MDV → upstream GUI contribution

SpatialData.js panel implementations are intentionally basic. MDV builds richer MUI
panels locally (reorderable accordions, histogram channel controls, table-aware shapes
panels). Once designs stabilise, contribute hardened headless primitives upstream rather
than opinionated styling.
