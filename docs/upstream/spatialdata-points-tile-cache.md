# spatialdata.js: the points tile cache holds 16 tiles, so zooming rebuilds them (draft)

Draft issue for `@spatialdata/core` / `@spatialdata/layers`, found in MDV's SpatialData chart.
Environment: Chrome (Dawn, macOS), `@spatialdata/*` 0.13.2 (the code is unchanged in 0.14.0),
`@deck.gl/core` 9.4.0, `@luma.gl/*` 9.4.2. Not posted.

## Summary

A tiled points layer caches far fewer tiles than one zoom touches. Zooming in or out over a
few levels therefore keeps evicting tiles and rebuilding them on the way back: a fresh
`getTileData`, a new `PointsScatterplotLayer` per tile, and with it a new luma `Model`
(on WebGPU, a full WGSL assembly — see `luma-wgsl-assembly-cache.md`). In MDV this showed
up as zooming dropping frames on a points-zarr dataset, even back and forth over the same
range.

`multiscaleTileGrid` (and `mortonTileGrid` the same way) sizes the cache as

```ts
const estimatedRowsPerTile = POINTS_TILE_TARGET_ROWS;           // 400_000
const budget = input.cacheRowBudget ?? DEFAULT_POINTS_MEMORY_CAP; // 4_000_000
const maxCacheSize = Math.min(512, Math.max(16, Math.round(budget / estimatedRowsPerTile)));
```

`4e6 / 4e5` is 10, clamped to 16. But `POINTS_TILE_TARGET_ROWS` is a ceiling the loader
aims below, not what a tile holds. On the dataset below, the cached tiles held 6 to 3,582
points each (most under 200), so the 4M-row budget the comment in `mortonTiledStrategy`
describes would cover hundreds of tiles; the cache is held to 16.

## Measurement

MDV project "mdv_atera_kidney": a `points-zarr-tiled` transcripts element, 1,186,813,834
rows, `finestTileSize` 100, deck grid `minZoom -5` / `maxZoom 2`. Sixty wheel steps
zooming about 4 levels, over a range already visited (so nothing new needs loading):

| Tile cache | Tile layers rebuilt | Frame median / p90 / max | Frames over 33ms |
|---|---|---|---|
| 16 (as computed) | 39–51 per sweep | 16.6 / 73.5 / 117ms | 11 of 60 |
| 512 (forced in the page) | 0 | 16.7 / 17.6 / 22ms | 0 of 60 |

77 tiles were resident after the sweeps with the larger cache. Each rebuild cost ~14ms of
main thread on WebGPU (layer init, Model, attribute upload), ~8ms on WebGL.

## Suggestion

Bound the cache by what tiles actually hold rather than by the per-tile target: deck's
`TileLayer` takes `maxCacheByteSize` and reads each tile's `content.byteLength`, so the
columnar tile content could report its byte size and the grid pass a byte budget derived
from `DEFAULT_POINTS_MEMORY_CAP`. That keeps the stated worst case ADR 0005 asks for while
letting sparse tiles stay. A smaller change: count rows per resident tile and evict against
the row budget, keeping `maxCacheSize` as a backstop at 512.

## Fix

The byte-bounded version is on SpatialData.ts branch `points-tile-cache-bytes` (33fc8dd,
not pushed): tiles report their batch `byteLength`, `maxCacheByteSize` is the row budget
at 16 bytes a row (64 MB), and `maxCacheSize` is a 128-tile backstop, kept low because
deck diffs a sublayer per cached tile every frame. MDV linked to that build, same
60-step sweeps on the dataset above, with the old bounds forced in the page for the
baseline:

| Bounds | Tiles rebuilt per warm sweep | Frame p90 | Frames over 33ms |
|---|---|---|---|
| 16 tiles, no byte bound | 55–61 | 45.6–46.9ms | 9–13 of 60 |
| 64 MB, 128-tile backstop | 0 | 17.3–17.7ms | 0 of 60 |
