import type { OrthographicViewState } from "@deck.gl/core";

export type DensityGridLayout = {
    columns: number;
    rows: number;
    /** Square cell edge length in pixels. */
    cellSize: number;
    cellWidth: number;
    cellHeight: number;
    contentWidth: number;
    contentHeight: number;
};

export type DensityGridCellBounds = {
    x: number;
    y: number;
    width: number;
    height: number;
};

export const DEFAULT_MIN_CELL_WIDTH = 260;

/** Chart types that expose the density grid setting (see `includeDensityModeToggle`). */
export const DENSITY_GRID_CHART_TYPES = new Set(["DeckContourScatter", "DeckDensity"]);

export function supportsDensityGridMode(chartType: string | undefined) {
    return typeof chartType === "string" && DENSITY_GRID_CHART_TYPES.has(chartType);
}

export function getDensityGridLayout(
    viewportWidth: number,
    _viewportHeight: number,
    fieldCount: number,
    cellSize = DEFAULT_MIN_CELL_WIDTH,
): DensityGridLayout {
    const safeFieldCount = Math.max(1, fieldCount);
    const safeViewportWidth = Math.max(1, viewportWidth);
    const safeCellSize = Math.max(1, cellSize);
    const preferredColumns = Math.max(1, Math.floor(safeViewportWidth / safeCellSize));
    const columns = Math.min(safeFieldCount, preferredColumns);
    const rows = Math.ceil(safeFieldCount / columns);

    return {
        columns,
        rows,
        cellSize: safeCellSize,
        cellWidth: safeCellSize,
        cellHeight: safeCellSize,
        contentWidth: columns * safeCellSize,
        contentHeight: rows * safeCellSize,
    };
}

export function getDensityGridCellBounds(layout: DensityGridLayout, index: number): DensityGridCellBounds {
    const row = Math.floor(index / layout.columns);
    const column = index % layout.columns;

    return {
        x: column * layout.cellWidth,
        y: row * layout.cellHeight,
        width: layout.cellWidth,
        height: layout.cellHeight,
    };
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
