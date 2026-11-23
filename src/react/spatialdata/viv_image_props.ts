import { findRenderStackEntry, type RenderStack } from "@spatialdata/layers";

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
