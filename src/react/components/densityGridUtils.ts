import type { Layer, OrthographicViewState } from "@deck.gl/core";
import {
    getChartArrayViewId,
    getChartArrayViewStates,
    hasUsableOrthographicViewState,
    getVivGridDetailViewId,
} from "./chartArrayGridUtils";
import { shouldDrawDeckLayerInViewport, type DeckLayerScopeInput } from "./deckLayerViewportScope";

export { getChartArrayViewStates, hasUsableOrthographicViewState, getVivGridDetailViewId };

/** Chart types that expose the density grid setting (see `includeDensityModeToggle`). */
export const DENSITY_GRID_CHART_TYPES = new Set([
    "DeckContourScatter",
    "DeckDensity",
    "VivMdvRegionReact",
]);

export function supportsDensityGridMode(chartType: string | undefined) {
    return typeof chartType === "string" && DENSITY_GRID_CHART_TYPES.has(chartType);
}

export function getDensityGridViewId(chartId: string, fieldId: string, index: number) {
    return getChartArrayViewId(chartId, fieldId, index, "density-grid");
}

export function matchesDensityGridView(layerId: string, viewId: string) {
    return layerId === viewId || layerId.startsWith(`${viewId}-`);
}

export function isDensityGridViewport(viewportId: string) {
    return viewportId.startsWith("density-grid-");
}

export function isEditableSelectionLayerId(layerId: string) {
    return layerId.startsWith("selection_");
}

export type DeckCanvasViewport = {
    id: string;
    x?: number;
    y?: number;
    width?: number;
    height?: number;
    unproject: (coords: number[]) => number[];
};

/** Viewport whose bounds contain the canvas point, or the topmost match when cells overlap in z-order. */
export function getViewportAtCanvasPoint(
    viewports: readonly DeckCanvasViewport[],
    canvasX: number,
    canvasY: number,
): DeckCanvasViewport | undefined {
    for (let i = viewports.length - 1; i >= 0; i--) {
        const vp = viewports[i];
        const x = vp.x ?? 0;
        const y = vp.y ?? 0;
        const width = vp.width ?? 0;
        const height = vp.height ?? 0;
        if (canvasX >= x && canvasX < x + width && canvasY >= y && canvasY < y + height) {
            return vp;
        }
    }
    return viewports[0];
}

/** Unproject canvas pixels through the sub-viewport under the pointer (multi-view grids). */
export function unprojectCanvasPoint(
    viewports: readonly DeckCanvasViewport[],
    canvasX: number,
    canvasY: number,
): number[] {
    const viewport = getViewportAtCanvasPoint(viewports, canvasX, canvasY);
    if (!viewport) return [0, 0, 0];
    const localX = canvasX - (viewport.x ?? 0);
    const localY = canvasY - (viewport.y ?? 0);
    return viewport.unproject([localX, localY]);
}

export type LayerViewportFilterInput = DeckLayerScopeInput;

/**
 * Deck layerFilter for MDVivViewer: single detail viewports use viv id suffixes;
 * density-grid cells route layers by {@link DeckLayerViewportScope} and optional viewId.
 */
export function shouldDrawLayerInViewport(
    layer: LayerViewportFilterInput,
    viewportId: string,
    vivIdForViewport: string,
) {
    return shouldDrawDeckLayerInViewport(layer, viewportId, {
        overlayDetailVivId: vivIdForViewport,
    });
}

/** layerFilter for Deck-only density grid (no viv id suffix on viewports). */
export function shouldDrawLayerInDeckDensityGrid(layer: LayerViewportFilterInput, viewportId: string) {
    return shouldDrawDeckLayerInViewport(layer, viewportId);
}

export function getDensityGridViewStates(
    viewIds: string[],
    sharedViewState: OrthographicViewState,
): Record<string, OrthographicViewState> {
    return getChartArrayViewStates(viewIds, sharedViewState);
}

type CloneableDeckLayer = Layer & {
    clone: (props: Record<string, unknown>) => Layer;
};

export function cloneDeckLayer(layer: CloneableDeckLayer, props: Record<string, unknown>): Layer {
    return layer.clone(props);
}

export function getSerializableViewState(viewState: OrthographicViewState): OrthographicViewState {
    return {
        target: viewState.target,
        zoom: viewState.zoom,
        minZoom: viewState.minZoom,
        maxZoom: viewState.maxZoom,
    };
}
