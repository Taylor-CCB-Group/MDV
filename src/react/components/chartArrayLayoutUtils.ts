/**
 * Layout helpers for {@link ChartArrayLayout} (prototype — heuristics and sizing may change).
 */

/** Layout mode for {@link ChartArrayLayout}; overridable later for user prefs. */
export type ChartArrayLayoutMode = "single" | "grid";

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
    // two cells go through the grid too, so a tall panel stacks them rather than splitting it in two
    return cellCount <= 1 ? "single" : "grid";
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
        parseCssPx(styles.getPropertyValue("--mdv-chart-array-min-cell")) || CHART_ARRAY_LAYOUT_DEFAULT_SIZING.minCell;

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

export type ChartArrayGridChoice = {
    columns: number;
    /** True when the cells fill the scrollport (no scrolling); otherwise they are square and scroll. */
    fillsViewport: boolean;
};

/**
 * Columns for a grid of `cellCount` cells in an `innerWidth` x `innerHeight` area (padding excluded).
 *
 * When every cell can be at least {@link ChartArrayLayoutSizing.minCell} in both directions without
 * scrolling, picks the column count that makes cells as large as possible in their *smaller*
 * dimension, so a wide panel gets a row and a tall one a column rather than a few tall columns.
 * Otherwise the cells are square and scroll, wrapping as {@link getPreferredGridColumnCount} does.
 */
export function chooseChartArrayGrid(
    cellCount: number,
    innerWidth: number,
    innerHeight: number,
    sizing: ChartArrayLayoutSizing = CHART_ARRAY_LAYOUT_DEFAULT_SIZING,
): ChartArrayGridChoice {
    if (cellCount <= 0 || innerWidth <= 0) return { columns: 1, fillsViewport: true };
    const maxColumns = Math.min(estimateGridColumnCount(innerWidth, sizing), cellCount);
    let best: { columns: number; side: number } | undefined;
    for (let columns = 1; columns <= maxColumns; columns++) {
        const rows = Math.ceil(cellCount / columns);
        const width = (innerWidth - (columns - 1) * sizing.gap) / columns;
        const height = (innerHeight - (rows - 1) * sizing.gap) / rows;
        const side = Math.min(width, height);
        if (side < sizing.minCell) continue;
        // more than half a pixel better; on a tie the fewer columns (and fewer empty slots) win
        if (!best || side > best.side + 0.5) best = { columns, side };
    }
    if (best) return { columns: best.columns, fillsViewport: true };
    return { columns: getPreferredGridColumnCount(cellCount, innerWidth, sizing), fillsViewport: false };
}
