# WebGPU picking: upstream issues (drafts)

Draft issues for deck.gl, SpatialData.js and deck.gl-community/editable-layers, found while
bringing MDV's SpatialData chart up on WebGPU. Each section is written to paste into an issue:
title, summary, standalone repro, what happens, the cause where known, and a suggested fix.

Every repro was run against a bare `Deck` (no MDV code). Environment:

- Chrome 152 (Dawn, `metal-3` adapter, macOS)
- `@deck.gl/core`, `@deck.gl/layers` 9.4.0; `@luma.gl/*` 9.4.2
- `@spatialdata/layers` 0.13.1
- `@deck.gl-community/editable-layers` 9.4.1; the 9.4.8 source was also checked

Each cause was also checked against deck.gl `master` at `b04c8f2bda` (2026-10-06), where it is
still present. No existing issue or PR reports any of these WebGPU failures. One is related:
deck.gl #9733 lists `pickMultipleObjectsAsync()` as a TODO (see issue 4).

| # | Repo | Issue |
|---|------|-------|
| 1 | deck.gl | WebGPU picks read the vertically mirrored pixel; point and rectangle |
| 2 | SpatialData.js | `applyWebGPUPickingFix` makes deck's rectangle picks throw |
| 3 | deck.gl | `PathLayer` WGSL has no picking output: picks return the line colour as an index |
| 4 | deck.gl | No multi-object (depth) pick on WebGPU |
| 5 | deck.gl-community | editable-layers picks synchronously; on WebGPU it throws and swallows other layers' clicks |

All repros share this setup:

```js
import { Deck, OrthographicView } from '@deck.gl/core';
import { ScatterplotLayer, PathLayer } from '@deck.gl/layers';
import { webgpuAdapter } from '@luma.gl/webgpu';

// 400×300 canvas; 1 world unit = 1 CSS pixel, world (0, 0) at the centre,
// so world (x, y) is drawn at CSS (x + 200, y + 150).
function makeDeck(layers, onDeviceInitialized) {
  return new Promise((resolve) => {
    const deck = new Deck({
      parent: document.body,
      width: 400,
      height: 300,
      views: new OrthographicView(),
      initialViewState: { target: [0, 0, 0], zoom: 0 },
      deviceProps: { type: 'webgpu', adapters: [webgpuAdapter] },
      onDeviceInitialized,
      layers,
      onAfterRender: () => setTimeout(() => resolve(deck), 300),
    });
  });
}

// A dot drawn at CSS (100, 50); its vertical mirror is (100, 250).
const dot = new ScatterplotLayer({
  id: 'dot',
  data: [{ position: [-100, -100] }],
  getPosition: (d) => d.position,
  getRadius: 15,
  radiusUnits: 'pixels',
  getFillColor: [255, 80, 0],
  pickable: true,
});
```

---

## 1. deck.gl: on WebGPU, picks read the vertically mirrored pixel (point and rectangle)

**Summary.** On a WebGPU device, `pickObjectAsync` and `pickObjectsAsync` find what is drawn at
the vertical mirror of the requested position. Hover, click and the public pick APIs are all
affected.

**Repro.**

```js
const deck = await makeDeck([dot]);

await deck.pickObjectAsync({ x: 100, y: 50, radius: 2 });   // null      (expected: the dot)
await deck.pickObjectAsync({ x: 100, y: 250, radius: 2 });  // the dot   (the mirror)

await deck.pickObjectsAsync({ x: 80, y: 30, width: 40, height: 40 });   // []
await deck.pickObjectsAsync({ x: 80, y: 230, width: 40, height: 40 });  // [the dot]
```

**Cause.** `DeckPicker` converts every pick position with
`canvasContext.cssToDevicePixels([x, y], true)`. The `true` asks for WebGL's bottom-left origin,
on WebGPU as well, where render targets are top-left. See `modules/core/src/lib/deck-picker.ts`
on master: lines 266 and 432 (point picks), and 588/595 and 694/701 (rectangle picks).

**Suggested fix.** Pass `device.type !== 'webgpu'` as the `yInvert` argument, as luma's own
picker does. The rectangle paths then also need their edges taken from the top-left:
`top = leftTop.y`, `bottom = rightBottom.y + rightBottom.height`, `height = bottom - top`. Today
they take `deviceTop = leftTop.y + leftTop.height` and `height = deviceTop - deviceBottom`,
which only holds with a bottom-left origin (see issue 2).

---

## 2. SpatialData.js: `applyWebGPUPickingFix` makes deck's rectangle picks throw

**Summary.** `applyWebGPUPickingFix` works around issue 1 by making `cssToDevicePixels` ignore
`yInvert`. That fixes point picks. But deck's rectangle-pick maths assumes a bottom-left origin,
so with the patch the rectangle gets a negative height, and every `pickObjectsAsync` throws.

**Repro.**

```js
import { applyWebGPUPickingFix } from '@spatialdata/layers';

const deck = await makeDeck([dot], (device) => applyWebGPUPickingFix(device));

await deck.pickObjectAsync({ x: 100, y: 50, radius: 2 });              // the dot (fixed)
await deck.pickObjectsAsync({ x: 80, y: 30, width: 40, height: 40 });  // throws
```

```
TypeError: Failed to execute 'createBuffer' on 'GPUDevice': Failed to read the 'size' property
from 'GPUBufferDescriptor': Value is outside the 'unsigned long long' value range.
    at WebGPUDevice.createBuffer
    at DeckPicker._readTextureDataAsync
    at DeckPicker._drawAndSampleAsync
    at DeckPicker._pickVisibleObjectsAsync
    at DeckPicker.pickObjectsAsync
    at Deck.pickObjectsAsync
```

**Cause.** `_pickVisibleObjectsAsync` computes
`deviceTop = leftTop.y + leftTop.height` and `deviceBottom = rightBottom.y`, so
`height = deviceTop - deviceBottom`. With the patch, `leftTop.y < rightBottom.y`, so the height
is negative.

**Suggested fix.** The real fix is upstream (issue 1). Until then, either document that
rectangle picks are unavailable while the patch is applied, or also patch the rectangle path
(`DeckPicker.prototype._pickVisibleObjectsAsync` and its synchronous twin) to compute a
top-left rectangle. MDV avoids rectangle picks for now (it uses point picks only).

---

## 3. deck.gl: `PathLayer`'s WGSL has no picking output; picks return the line colour as an index

**Summary.** On WebGPU, a pick on a `PathLayer` returns `picked: true` with no `object`, and an
`index` that is the line's display colour decoded as a picking colour. This also affects layers
built from `PathLayer`: `GeoJsonLayer` and `PolygonLayer` outlines (confirmed on
`GeoJsonLayer`'s `polygons-stroke`), and editable-layers' `EditablePathLayer`.

**Repro.** Pick at the mirrored position because of issue 1, or apply the issue 2 workaround
and pick at (200, 50); the result is the same either way.

```js
const path = new PathLayer({
  id: 'path',
  data: [{ path: [[-150, -100], [150, -100]] }],  // horizontal, at CSS y = 50
  getPath: (d) => d.path,
  getWidth: 10,
  widthUnits: 'pixels',
  getColor: [10, 20, 30],
  pickable: true,
});
const deck = await makeDeck([path]);

const info = await deck.pickObjectAsync({ x: 200, y: 250, radius: 2 });
// { picked: true, index: 1971209, object: undefined, ... }
// 1971209 + 1 = 0x1E140A = rgb(10, 20, 30): the line colour, not a picking colour.
```

With `getColor: [150, 150, 150]` the index is 9868949 (`0x969696`), so the index follows the
colour.

**Cause.** In `modules/layers/src/path-layer/path-layer.wgsl.ts`, the vertex shader computes
`geometry.pickingColor` (line 238 on master), but `fragmentMain` never outputs it. There is no
`picking.isActive` branch, and `DECKGL_FILTER_COLOR` is a comment in WGSL. So the picking pass
draws the path in its display colour. Compare `scatterplot-layer.wgsl.ts`, which handles
picking explicitly:

```wgsl
if (picking.isActive > 0.5) {
  if (!picking_isColorValid(varyings.pickingColor)) {
    discard;
  }
  return vec4<f32>(varyings.pickingColor, 1.0);
}
```

**Suggested fix.** Carry the picking colour as a varying and add the same branch to
`PathLayer`'s fragment shader. Do it after the analytic-antialiasing discards and before
`deckgl_premultiplied_alpha`. Of the 9.4 WGSL layers, `column-layer.wgsl.ts` is the only other
one without a picking branch; it is probably affected too, but this was not tested at runtime.

---

## 4. deck.gl: no multi-object (depth) pick on WebGPU

**Summary.** On WebGPU there is no way to get every object under a point. `pickMultipleObjects`
is documented as WebGL-only, and on WebGPU it throws. There is no `pickMultipleObjectsAsync`.
`pickObjectsAsync` returns only the top-most object at each pixel of a rectangle.

**Repro.**

```js
const deck = await makeDeck([dot]);
deck.pickMultipleObjects({ x: 100, y: 250 });
// Error: not implemented   (WebGPUDevice.readPixelsToArrayWebGL, from DeckPicker._drawAndSample)
```

**Why it matters.**

- Multi-layer tooltips: everything under the cursor (image channels, shapes, points, gates).
- Editing layers: an edit handle and the feature beneath it (issue 5).
- SpatialData.js's aggregate hover currently falls back to a single pick on WebGPU for this
  reason.

**Request.** `pickMultipleObjectsAsync({ x, y, radius, depth, layerIds })`, the async counterpart
of `pickMultipleObjects`, as anticipated in the TODO list of #9733 ("Determine how to correctly
implement `pickMultipleObjectsAsync()`"). The open question there, whether it must snapshot
deck props between its successive picking passes, could be answered by accepting results
against the layer state at the time of the call.

---

## 5. deck.gl-community/editable-layers: synchronous picking throws on WebGPU and swallows other layers' clicks

**Summary.** `EditableLayer.getPicks` calls `deck.pickMultipleObjects`, which throws on WebGPU
(issue 4). It is called from the layer's pointer handlers on every `pointermove`, `click` and
pan event. The throw escapes the mjolnir handler, and a click is not dispatched to other layers
at all: a pickable layer's `onClick` stops firing as soon as an `EditableGeoJsonLayer` is
present.

**Repro.**

```js
import { EditableGeoJsonLayer, TranslateMode } from '@deck.gl-community/editable-layers';

let clicks = 0;
const clickableDot = new ScatterplotLayer({
  id: 'dot',
  data: [{ position: [-100, -100] }],
  getPosition: (d) => d.position,
  getRadius: 15,
  radiusUnits: 'pixels',
  getFillColor: [255, 80, 0],
  pickable: true,
  onClick: () => { clicks++; return true; },
});
const editable = new EditableGeoJsonLayer({
  id: 'edit',
  data: {
    type: 'FeatureCollection',
    // a rectangle away from the dot
    features: [{ type: 'Feature', properties: {}, geometry: { type: 'Polygon',
      coordinates: [[[50, 50], [150, 50], [150, 120], [50, 120], [50, 50]]] } }],
  },
  mode: TranslateMode,
  selectedFeatureIndexes: [0],
  pickable: true,
});

// The workaround from issue 2, so the dot is picked where it is drawn.
const deck = await makeDeck([clickableDot, editable], (device) => applyWebGPUPickingFix(device));

// Click the dot: pointermove, pointerdown, pointerup at CSS (100, 50) on the canvas.
```

| Layers | `onClick` calls | Uncaught errors |
|---|---|---|
| dot only | 1 | none |
| dot + `EditableGeoJsonLayer` | 0 | `Error: not implemented` |

`dist/editable-layers/editable-layer.js` still calls `pickMultipleObjects` in 9.4.8 (line 397).

**Suggested fix.**

1. Never let a pick failure escape the event handler: catch it, and treat it as no picks.
2. On WebGPU, do not rely on GPU picks taken at the event. An async pick cache (pick on
   `pointermove`, use the latest result when a drag starts) answers for where the pointer
   *was*. In practice, translate-versus-draw became a race when the pointer moved just before
   the press. MDV tried this first.
3. What MDV ships instead is a CPU pick for the editable layer's own features and edit handles:
   point-in-polygon and edge distance for features, distance for handles, all in screen space.
   An editing layer holds little geometry, so this is cheap, synchronous and exact. It also
   sidesteps issues 2 and 3, which bit the GPU-pick attempt (rectangle picks threw, and outline
   picks returned garbage indices with no `object`, which crashed `ModifyMode.getGuides` on
   `pick.object.geometry`). Something like this inside editable-layers, used when the device is
   WebGPU, would remove the need for host workarounds.

---

## MDV workarounds that these would retire

| Workaround | Where | Retired by |
|---|---|---|
| CPU picks for the selection layer on WebGPU | `src/editable-layers/cpuPicks.ts`, `MonkeyPatchEditableGeoJsonLayer.getPicks` in `src/lib/deckMonkeypatch.ts` | 5 (and 3, 4) |
| The selection shape is always selected, not selected on hover (hover selection raced async picks) | `selectedFeatureIndexes` in `src/react/spatial_context.tsx` | 5, or a multi-shape selection model |
| SpatialData.js's `applyWebGPUPickingFix` (applied by its viewers) | `@spatialdata/layers` | 1 |
