# Layer Group Portability Notes (MDV -> SpatialData.js)

This note captures likely extraction friction for moving the layer-group implementation from MDV into `SpatialData.js`.

## Can likely copy as-is

- `src/webgl/layerGroups.ts` contracts (`LayerGroup`, `LayerGroupScope`, composer).
- Cache key strategy from `src/webgl/chartArrayRenderCache.ts` (`getChartArrayStaticCacheKey`).
- Contour backend interface from `src/webgl/contourBackend.ts`.

## Needs adapter boundaries in MDV

- React/MobX config and hooks (`DeckDensityGridComponent`, chart config mutation) are MDV-specific.
- Viewport ID conventions (`density-grid-*`, viv suffixes) are app-specific and should remain adapter logic.
- Selection/gating layer IDs and scope tagging should stay outside reusable package interfaces.

## Deck/Luma compatibility checks

- MDV currently uses deck.gl `9.1.x` + luma.gl `9.1.x`; verify identical major/minor expectations in `SpatialData.js` before extraction.
- Ensure `device.createTexture`/`device.createFramebuffer` usage remains valid in target runtime.
- If API deltas are found during extraction, record them in a dedicated compatibility note before proceeding.

## Suggested package-level interfaces before porting

- `LayerGroup` and `LayerGroupComposer` as pure deck-agnostic orchestration types.
- `ContourBackend` as pluggable rendering backend contract.
- `StaticPassCache` as lifecycle-managed resource object with explicit `dispose()`.
