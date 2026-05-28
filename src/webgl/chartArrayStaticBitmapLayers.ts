import { BitmapLayer } from "@deck.gl/layers";
import { OrthographicViewport, type Layer, type OrthographicViewState } from "@deck.gl/core";
import type { Texture } from "@luma.gl/core";
import { tagDeckLayerViewportScope } from "@/react/components/deckLayerViewportScope";

/** Match Viv `DetailView` / deck `OrthographicView` defaults (not `flipY: false`). */
export const CHART_ARRAY_STATIC_VIEW_FLIP_Y = true;

/** World-space bounds for a full-canvas orthographic view (deck flipY=false). */
export function getOrthographicWorldBounds(
    viewState: OrthographicViewState,
    width: number,
    height: number,
): [number, number, number, number] {
    const target = viewState.target ?? [0, 0, 0];
    const zoom = resolveOrthographicZoom(viewState.zoom);
    const scale = 2 ** zoom;
    const worldWidth = width / scale;
    const worldHeight = height / scale;
    const left = target[0] - worldWidth / 2;
    const right = target[0] + worldWidth / 2;
    const bottom = target[1] - worldHeight / 2;
    const top = target[1] + worldHeight / 2;
    return [left, bottom, right, top];
}

function resolveOrthographicZoom(zoom: OrthographicViewState["zoom"]): number {
    if (Array.isArray(zoom)) {
        return Math.min(zoom[0], zoom[1]);
    }
    return Number.isFinite(zoom) ? Number(zoom) : 0;
}

/**
 * Corner bounds for BitmapLayer when sampling an FBO color attachment.
 * Uses the same {@link OrthographicViewport} math as the live grid viewports.
 */
export function getFramebufferCompositeBounds(
    viewState: OrthographicViewState,
    width: number,
    height: number,
): [[number, number], [number, number], [number, number], [number, number]] {
    const viewport = new OrthographicViewport({
        ...viewState,
        width,
        height,
        flipY: CHART_ARRAY_STATIC_VIEW_FLIP_Y,
        id: "chart-array-static-bounds",
    });
    const topLeft = viewport.unproject([0, 0]);
    const topRight = viewport.unproject([width, 0]);
    const bottomLeft = viewport.unproject([0, height]);
    const bottomRight = viewport.unproject([width, height]);
    // FBO color attachment uses WebGL origin; BitmapLayer shaders assume image-style V.
    return [
        [topLeft[0], topLeft[1]],
        [bottomLeft[0], bottomLeft[1]],
        [bottomRight[0], bottomRight[1]],
        [topRight[0], topRight[1]],
    ];
}

/** Composite the static shared-geometry texture into one grid viewport. */
export function createChartArrayStaticCompositeLayer(
    id: string,
    viewId: string,
    texture: Texture,
    bounds:
        | [number, number, number, number]
        | [[number, number], [number, number], [number, number], [number, number]],
): Layer {
    return tagDeckLayerViewportScope(
        new BitmapLayer({
            id,
            image: texture,
            bounds,
            alphaCutoff: 0,
            pickable: false,
        }),
        "per-viewport",
        { viewId },
    );
}
