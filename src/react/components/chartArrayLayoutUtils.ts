/**
 * Layout helpers for {@link ChartArrayLayout} (prototype — heuristics and sizing may change).
 */

/** Layout mode for {@link ChartArrayLayout}; overridable later for user prefs. */
export type ChartArrayLayoutMode = "single" | "pair" | "grid";

/** Mirrors defaults in chartArrayLayout.css. */
export type ChartArrayLayoutSizing = {
    gap: number;
    paddingX: number;
    paddingY: number;
    minCell: number;
};

export const CHART_ARRAY_LAYOUT_DEFAULT_SIZING: ChartArrayLayoutSizing = {
    gap: 8,
    paddingX: 24,
    paddingY: 24,
    minCell: 260,
};

export function getChartArrayLayoutMode(cellCount: number): ChartArrayLayoutMode {
    if (cellCount <= 1) return "single";
    if (cellCount === 2) return "pair";
    return "grid";
}

function parseCssPx(value: string): number {
    const parsed = Number.parseFloat(value);
    return Number.isFinite(parsed) ? parsed : 0;
}

export function readChartArraySizingFromRoot(root: HTMLElement): ChartArrayLayoutSizing {
    const styles = getComputedStyle(root);
    const paddingLeft = parseCssPx(styles.paddingLeft);
    const paddingRight = parseCssPx(styles.paddingRight);
    const paddingTop = parseCssPx(styles.paddingTop);
    const paddingBottom = parseCssPx(styles.paddingBottom);
    const gap =
        parseCssPx(styles.columnGap) ||
        parseCssPx(styles.rowGap) ||
        parseCssPx(styles.gap) ||
        CHART_ARRAY_LAYOUT_DEFAULT_SIZING.gap;
    const minCell =
        parseCssPx(styles.getPropertyValue("--mdv-chart-array-min-cell")) ||
        CHART_ARRAY_LAYOUT_DEFAULT_SIZING.minCell;

    return {
        gap,
        paddingX: paddingLeft + paddingRight,
        paddingY: paddingTop + paddingBottom,
        minCell,
    };
}

/** Matches `repeat(auto-fill, minmax(min(minCell, 100%), 1fr))` column count. */
export function estimateGridColumnCount(
    innerWidth: number,
    sizing: ChartArrayLayoutSizing = CHART_ARRAY_LAYOUT_DEFAULT_SIZING,
): number {
    if (innerWidth <= 0) return 1;
    const trackMin = Math.min(sizing.minCell, innerWidth);
    return Math.max(1, Math.floor((innerWidth + sizing.gap) / (trackMin + sizing.gap)));
}

/**
 * Column count that lets items share row width evenly, wrapping only when a single-row
 * share would be narrower than {@link ChartArrayLayoutSizing.minCell}.
 */
export function getPreferredGridColumnCount(
    cellCount: number,
    innerWidth: number,
    sizing: ChartArrayLayoutSizing = CHART_ARRAY_LAYOUT_DEFAULT_SIZING,
): number {
    if (cellCount <= 0) return 1;
    const maxColumns = Math.min(estimateGridColumnCount(innerWidth, sizing), cellCount);
    const equalShareWidth = (innerWidth - (cellCount - 1) * sizing.gap) / cellCount;
    if (equalShareWidth >= sizing.minCell) {
        return cellCount;
    }
    const rowsIfMaxColumns = Math.ceil(cellCount / maxColumns);
    return Math.max(1, Math.ceil(cellCount / rowsIfMaxColumns));
}

export function estimateSquareGridContentHeight(
    cellCount: number,
    innerWidth: number,
    sizing: ChartArrayLayoutSizing = CHART_ARRAY_LAYOUT_DEFAULT_SIZING,
): number {
    if (cellCount <= 0) return 0;
    const columns = getPreferredGridColumnCount(cellCount, innerWidth, sizing);
    const rows = Math.ceil(cellCount / columns);
    const cellWidth = (innerWidth - (columns - 1) * sizing.gap) / columns;
    return rows * cellWidth + (rows - 1) * sizing.gap + sizing.paddingY;
}

/** True when square cells at the current width would fit inside the scrollport without scrolling. */
export function chartArrayGridFitsViewport(
    cellCount: number,
    viewportHeight: number,
    innerWidth: number,
    sizing: ChartArrayLayoutSizing = CHART_ARRAY_LAYOUT_DEFAULT_SIZING,
): boolean {
    if (cellCount <= 0 || viewportHeight <= 0 || innerWidth <= 0) return true;
    return estimateSquareGridContentHeight(cellCount, innerWidth, sizing) <= viewportHeight;
}
