import type { Layer, OrthographicViewState } from "@deck.gl/core";
import { getConcreteFieldNames } from "@/charts/chartConfigUtils";
import type { FieldSpecs } from "@/lib/columnTypeHelpers";
import {
    getChartArrayViewId,
    getChartArrayViewStates,
    hasUsableOrthographicViewState,
    getVivGridDetailViewId,
    supportsChartArrayGridMode,
} from "./chartArrayGridUtils";
import { shouldDrawDeckLayerInViewport, type DeckLayerScopeInput } from "./deckLayerViewportScope";

export { getChartArrayViewStates, hasUsableOrthographicViewState, getVivGridDetailViewId };

/** @deprecated Use {@link supportsChartArrayGridMode} */
export const DENSITY_GRID_CHART_TYPES = new Set([
    "DeckContourScatter",
    "DeckDensity",
    "VivMdvRegionReact",
]);

/** @deprecated Use {@link supportsChartArrayGridMode} */
export function supportsDensityGridMode(chartType: string | undefined) {
    return supportsChartArrayGridMode(chartType);
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
    isPicking?: boolean,
) {
    return shouldDrawDeckLayerInViewport(layer, viewportId, {
        overlayDetailVivId: vivIdForViewport,
        isPicking,
    });
}

/** layerFilter for Deck-only density grid (no viv id suffix on viewports). */
export function shouldDrawLayerInDeckDensityGrid(
    layer: LayerViewportFilterInput,
    viewportId: string,
    isPicking?: boolean,
) {
    return shouldDrawDeckLayerInViewport(layer, viewportId, { isPicking });
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

export const CHART_ARRAY_OVERLAY_DETAIL_SUFFIX = "detail-react";
export const CHART_ARRAY_GRID_LAYER_SUFFIX = "chart-array-grid";

/** Overlay detail layers use `-#<chartId>detail-react#`; grid clones use `chart-array-grid` instead. */
export function toChartArrayGridLayerId(layerId: string): string {
    if (layerId.includes(CHART_ARRAY_GRID_LAYER_SUFFIX)) {
        return layerId;
    }
    if (layerId.includes(CHART_ARRAY_OVERLAY_DETAIL_SUFFIX)) {
        return layerId.replace(CHART_ARRAY_OVERLAY_DETAIL_SUFFIX, CHART_ARRAY_GRID_LAYER_SUFFIX);
    }
    return `${layerId}-${CHART_ARRAY_GRID_LAYER_SUFFIX}`;
}

export function cloneLayerForChartArrayGrid(
    layer: CloneableDeckLayer,
    props: Record<string, unknown> = {},
): Layer {
    return cloneDeckLayer(layer, {
        ...props,
        id: toChartArrayGridLayerId(layer.id),
    });
}

export function cloneDeckLayerForRender(
    layer: CloneableDeckLayer,
    props: Record<string, unknown> = {},
): Layer {
    return cloneDeckLayer(layer, {
        ...props,
        id: layer.id,
    });
}

export function cloneDeckLayer(layer: CloneableDeckLayer, props: Record<string, unknown>): Layer {
    return layer.clone(props);
}

export const DENSITY_GRID_EMPTY_STATE_MESSAGES = {
    noCellsConfigured: "Choose density fields to build the grid.",
    loadingCells: "Loading density fields...",
    noRows: "No rows remain after the current filters.",
} as const;

/** Count configured density columns (not loaded field count). */
export function getConfiguredDensityFieldCount(densityFields: FieldSpecs | undefined): number {
    if (densityFields === undefined || densityFields === null) return 0;
    return getConcreteFieldNames(densityFields).length;
}

export function getSerializableViewState(viewState: OrthographicViewState): OrthographicViewState {
    return {
        target: viewState.target,
        zoom: viewState.zoom,
        minZoom: viewState.minZoom,
        maxZoom: viewState.maxZoom,
    };
}

/** Normalize DeckGL view-state updates from single- or multi-view controllers. */
export function applyDeckViewStateChange(
    update: OrthographicViewState | Record<string, OrthographicViewState>,
    current?: OrthographicViewState,
): OrthographicViewState {
    if (
        update &&
        typeof update === "object" &&
        "target" in update &&
        Array.isArray((update as OrthographicViewState).target)
    ) {
        return getSerializableViewState(update as OrthographicViewState);
    }
    const record = update as Record<string, OrthographicViewState>;
    const next = Object.values(record).find(
        (viewState) =>
            viewState &&
            Array.isArray(viewState.target) &&
            Number.isFinite(Number(viewState.zoom)),
    );
    if (next) return getSerializableViewState(next);
    if (current) return getSerializableViewState(current);
    return { target: [0, 0, 0], zoom: 0 };
}
