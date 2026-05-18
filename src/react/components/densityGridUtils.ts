import type { Layer, OrthographicViewState } from "@deck.gl/core";
import {
    getChartArrayViewId,
    getChartArrayViewStates,
    hasUsableOrthographicViewState,
    getVivGridDetailViewId,
} from "./chartArrayGridUtils";

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

type LayerViewportFilterInput = {
    id: string;
    props?: { viewId?: string };
};

/**
 * Deck layerFilter for MDVivViewer: single detail viewports use viv id suffixes;
 * density-grid cells route layers by viewId / grid id prefix instead.
 */
export function shouldDrawLayerInViewport(
    layer: LayerViewportFilterInput,
    viewportId: string,
    vivIdForViewport: string,
) {
    if (layer.id.includes(vivIdForViewport)) {
        return true;
    }
    if (!isDensityGridViewport(viewportId)) {
        return false;
    }
    const layerViewId = layer.props?.viewId;
    if (typeof layerViewId === "string") {
        return (
            matchesDensityGridView(layer.id, layerViewId) ||
            matchesDensityGridView(layer.id, viewportId) ||
            layerViewId === viewportId
        );
    }
    return matchesDensityGridView(layer.id, viewportId);
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
