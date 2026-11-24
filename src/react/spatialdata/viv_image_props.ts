import {
    findRenderStackEntry,
    type RenderStack,
    type ResolvedAlphaFromIntensity,
    resolveAlphaFromIntensity,
} from "@spatialdata/layers";

function isRecord(value: unknown): value is Record<string, unknown> {
    return typeof value === "object" && value !== null && !Array.isArray(value);
}

/**
 * The Viv props an image entry is drawn with: the ones saved on its render stack entry, and
 * `pickable: false`.
 *
 * Images are drawn, never picked; no feature comes of an image pick. Viv's image layers are
 * pickable unless told otherwise, and SpatialData.js's `DeviceAdaptiveImageLayer` draws its
 * Viv layer with the Viv layer's own `pickable`, not its own (which is off). On WebGPU the
 * image tiles then fill the picking buffer wherever no point covers a pixel: a click a pixel
 * off a point picked the image, which has no click handler, so the row was not highlighted.
 * WebGL's Viv tiles write no picking colour, which hid it there.
 */
export function vivImageLayerProps(stack: RenderStack | undefined, layerId: string): Record<string, unknown> {
    const entry = stack ? findRenderStackEntry(stack, layerId) : undefined;
    const saved = entry?.kind === "spatial" ? entry.props.vivLayerProps : undefined;
    return { ...(isRecord(saved) ? saved : {}), pickable: false };
}

/**
 * An image entry's "opacity from brightness" (`props.alphaFromIntensity`), or null when it is
 * off. On WebGL it needs `AlphaFromIntensityExtension` on that layer alone: the extension is
 * not neutral when the option is off (it falls back to threshold 0, gain 1, which discards
 * black and rescales colour), so it must not be on every image layer. WebGPU's raster layer
 * reads the same prop itself.
 */
export function imageAlphaFromIntensity(
    stack: RenderStack | undefined,
    layerId: string,
): ResolvedAlphaFromIntensity | null {
    const entry = stack ? findRenderStackEntry(stack, layerId) : undefined;
    return entry?.kind === "spatial" ? resolveAlphaFromIntensity(entry.props.alphaFromIntensity) : null;
}
