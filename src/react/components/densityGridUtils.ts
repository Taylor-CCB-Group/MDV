import type { OrthographicViewState } from "@deck.gl/core";

/** Chart types that expose the density grid setting (see `includeDensityModeToggle`). */
export const DENSITY_GRID_CHART_TYPES = new Set(["DeckContourScatter", "DeckDensity"]);

export function supportsDensityGridMode(chartType: string | undefined) {
    return typeof chartType === "string" && DENSITY_GRID_CHART_TYPES.has(chartType);
}

export function getDensityGridViewId(chartId: string, fieldId: string, index: number) {
    const safeChartId = chartId.replace(/[^A-Za-z0-9_-]/g, "_");
    const safeFieldId = fieldId.replace(/[^A-Za-z0-9_-]/g, "_");
    return `density-grid-${safeChartId}-${index}-${safeFieldId}`;
}

export function matchesDensityGridView(layerId: string, viewId: string) {
    return layerId === viewId || layerId.startsWith(`${viewId}-`);
}

export function getDensityGridViewStates(
    viewIds: string[],
    sharedViewState: OrthographicViewState,
): Record<string, OrthographicViewState> {
    return Object.fromEntries(
        viewIds.map((viewId) => [
            viewId,
            {
                ...sharedViewState,
            },
        ]),
    );
}

export function hasUsableOrthographicViewState(viewState: OrthographicViewState | undefined): boolean {
    if (!viewState?.target || viewState.target.length < 2) return false;
    const x = Number(viewState.target[0]);
    const y = Number(viewState.target[1]);
    if (!Number.isFinite(x) || !Number.isFinite(y)) return false;
    return Number.isFinite(Number(viewState.zoom));
}

