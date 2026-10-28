# intraspatial: DensityLayer destroys buffers a splat in flight still uses (draft)

Draft issue for `@intraspatial/deck` / `intraspatial`, found in MDV's SpatialData chart on
WebGPU. Environment: Chrome (Dawn, `metal-3`, macOS), `@intraspatial/deck` 0.1.0,
`@deck.gl/core` 9.4.0, `@luma.gl/*` 9.4.2.

## Summary

When deck finalizes a `DensityLayer` while its splat is still in flight, the layer destroys
the raw `GPUBuffer`s that splat goes on to use. Each one produces an uncaptured validation
error naming an unlabeled buffer:

```
[Buffer (unlabeled)] used in submit while destroyed.
 - While calling [Queue].WriteBuffer([Buffer (unlabeled)], (0 bytes), data, (4 bytes))
[Buffer (unlabeled)] used in submit while destroyed.
 - While calling [Queue].Submit([[CommandBuffer]])
```

It happens whenever a host replaces density layers mid-splat. MDV does this when its contour
field list changes, and it has up to 60 `DensityLayer`s per chart sharing one queue, which
widens each splat's window. Because the buffers are unlabeled, the host cannot tell which
layer the error belongs to.

## Repro

A bare `Deck` on WebGPU, 40 layers of 300k points, replaced by layers with new ids every
40 ms. The page has to be in the foreground (deck needs `requestAnimationFrame` to process
layer changes). `Probe` only records whether the layer was busy when finalized.

```js
import { Deck, OrthographicView } from "@deck.gl/core";
import { webgpuAdapter } from "@luma.gl/webgpu";
import { DensityLayer } from "@intraspatial/deck";

const log = [];
class Probe extends DensityLayer {
    finalizeState(ctx) {
        log.push(Boolean(this.state?.busy));
        super.finalizeState(ctx);
    }
}
Probe.layerName = "Probe";

const errors = [];
const n = 300_000;
const x = Float32Array.from({ length: n }, Math.random);
const y = Float32Array.from({ length: n }, Math.random);
let ready;
const init = new Promise((r) => (ready = r));
const deck = new Deck({
    parent: document.body,
    deviceProps: { type: "webgpu", adapters: [webgpuAdapter], onError: (e) => (errors.push(e.message), true) },
    views: new OrthographicView(),
    initialViewState: { target: [0.5, 0.5, 0], zoom: 9 },
    layers: [],
    onDeviceInitialized: ready,
});
await init;
for (let gen = 0; gen < 6; gen++) {
    deck.setProps({
        layers: Array.from({ length: 40 }, (_, k) => new Probe({ id: `g${gen}-${k}`, x, y, bandwidth: 0.01 + k * 0.0005 })),
    });
    await new Promise((r) => setTimeout(r, 40));
}
deck.setProps({ layers: [] });
await new Promise((r) => setTimeout(r, 4000));
console.log({ finalized: log.length, finalizedWhileBusy: log.filter(Boolean).length, errors: errors.length });
```

Observed: `{ finalized: 240, finalizedWhileBusy: 40, errors: 200 }`, all the two messages
above. With one small layer the splat finishes within a frame and nothing is reported, so a
single-layer test passes.

## Cause

`finalizeState` destroys `field`, `max` and `points` unconditionally. `splat()` keeps running
after the layer is finalized: after `await splatDensityResident(...)` it calls
`fieldMaxResident(s.field, cells, s.max, ...)`, which, after its own `await getCtx(device)`,
does `queue.writeBuffer(dst, 0, new Uint32Array([0]))`: the 4-byte write at offset 0 in the
first message. `splatDensityResident` has internal awaits before its submit too, so `points`
and `field` can be destroyed before that submit as well.

`uploadPoints` already handles this case for `points` (it retires the buffer while `busy`);
`finalizeState` doesn't.

## Suggested fix

Make it a rule that a raw buffer is destroyed only while no splat holds it, so the error
cannot happen rather than being unlikely:

- In `finalizeState`, set `s.finalized = true`. If `s.busy`, move `field`, `max` and
  `points` to `s.retired` instead of destroying them; the existing `finally` in `refresh()`
  destroys them when the splat settles. (The luma wrappers and the model can still go
  immediately: luma never destroys a handle it was given.)
- In `splat()`, return after each `await` if `s.finalized`, so a finalized layer submits
  no further work. In `refresh()`'s `finally`, skip `again`/`setNeedsRedraw` once finalized.

A test can then pin it: finalize while `busy` (the repro above, or a stubbed
`splatDensityResident` that resolves on demand) and assert no device errors.

## Labels

Give every raw resource a `label`, so validation errors name their owner:

- `DensityLayer`: `${this.props.id}:field`, `:max`, `:points`, `:ramp`. With the layer id in
  the label, hosts can attribute an error to the layer that caused it. MDV's
  `entryForLabel` in `src/react/spatialdata/layer_failures.ts` already does this for luma
  resources, whose labels are their ids.
- The kernels on this path: `fieldMax`'s uniform buffer and `splatDensity`'s `ptsBuf`,
  `rbRaw` and `uniBuf`, named after the kernel. Many other intraspatial kernels create
  unlabeled buffers and textures too; a lint rule or a `createBuffer` wrapper that requires
  a label would keep it that way.
