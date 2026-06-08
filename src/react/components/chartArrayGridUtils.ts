import type { OrthographicViewState } from "@deck.gl/core";

const CHART_ARRAY_GRID_CHART_TYPES = new Set([
    "DeckContourScatter",
    "DeckDensity",
    "VivMdvRegionReact",
]);

export function getChartArrayViewId(
    chartId: string,
    cellKey: string,
    index: number,
    prefix = "chart-array-grid",
) {
    const safeChartId = chartId.replace(/[^A-Za-z0-9_-]/g, "_");
    const safeCellKey = cellKey.replace(/[^A-Za-z0-9_-]/g, "_");
    return `${prefix}-${safeChartId}-${index}-${safeCellKey}`;
}

export function getChartArrayViewStates(
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

export function getVivGridDetailViewId(chartId: string, cellKey: string, index: number) {
    return `${getChartArrayViewId(chartId, cellKey, index, "density-grid")}detail-react`;
}

export function supportsChartArrayGridMode(chartType: string | undefined): boolean {
    return typeof chartType === "string" && CHART_ARRAY_GRID_CHART_TYPES.has(chartType);
}

export function isChartArrayGridMode({
    chartType,
    dimension,
    layoutMode,
    cellCount,
}: {
    chartType: string | undefined;
    dimension: string | undefined;
    layoutMode: string | undefined;
    cellCount: number;
}): boolean {
    return (
        supportsChartArrayGridMode(chartType) &&
        dimension === "2d" &&
        layoutMode === "grid" &&
        cellCount > 0
    );
}
