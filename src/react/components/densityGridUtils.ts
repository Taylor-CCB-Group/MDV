import type { OrthographicViewState } from "@deck.gl/core";

export type DensityGridLayout = {
    columns: number;
    rows: number;
    cellWidth: number;
    cellHeight: number;
};

export type DensityGridCellBounds = {
    x: number;
    y: number;
    width: number;
    height: number;
};

const DEFAULT_MIN_CELL_WIDTH = 260;

export function getDensityGridLayout(
    width: number,
    height: number,
    fieldCount: number,
    minCellWidth = DEFAULT_MIN_CELL_WIDTH,
): DensityGridLayout {
    const safeFieldCount = Math.max(1, fieldCount);
    const safeWidth = Math.max(1, width);
    const safeHeight = Math.max(1, height);
    const preferredColumns = Math.max(1, Math.floor(safeWidth / Math.max(1, minCellWidth)));
    const columns = Math.min(safeFieldCount, preferredColumns);
    const rows = Math.ceil(safeFieldCount / columns);

    return {
        columns,
        rows,
        cellWidth: safeWidth / columns,
        cellHeight: safeHeight / rows,
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
