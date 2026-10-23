import { type MeasureText, type XLabelLayout, layoutFixedXLabels, layoutXLabels } from "@/charts/dotplot/xAxisLabels";
import type { DotPlotXAxisConfig } from "./dotPlotConfig";
import type { DotPlotGroup } from "./dotPlotStats";

export type RowOrder = "data" | "alphabetical";

/**
 * X labels for one panel, which suits them all (same width, same fields), as the legacy dot
 * plot lays them out: fitted with `fit_labels`, otherwise at the axis settings, shortened or
 * thinned to suit them.
 */
export function getXLabelLayout(
    axis: DotPlotXAxisConfig,
    labels: string[],
    space: { plotLeft: number; plotWidth: number; chartHeight: number },
    measure: MeasureText,
): XLabelLayout {
    const angle = axis.rotate_labels ? 45 : 0;
    if (labels.length === 0 || space.plotWidth <= 0) {
        return { angle, fontSize: axis.tickfont, margin: axis.size, labels };
    }
    const { tickfont: fontSize, truncate_labels: truncate, thin_labels: thin } = axis;
    const input = { labels, ...space, fontSize, titleSize: 0, truncate, thin };
    return axis.fit_labels
        ? layoutXLabels(input, measure)
        : layoutFixedXLabels({ ...input, angle, margin: axis.size }, measure);
}

/** Space above the plot area, as in the legacy dot plot. */
export const DOT_PLOT_TOP_MARGIN = 10;
export const DOT_PLOT_RIGHT_MARGIN = 10;
/** Largest dot radius in pixels, as in the legacy dot plot. */
export const DOT_PLOT_MAX_RADIUS = 25;

/**
 * Row category codes to draw, in display order. A row is drawn if it has data in at least one
 * panel, so every panel shows the same rows.
 */
export function getVisibleRows(panels: DotPlotGroup[][], values: string[], order: RowOrder): number[] {
    const rows = values.map((_, r) => r).filter((r) => panels.some((groups) => (groups[r]?.count ?? 0) > 0));
    if (order === "alphabetical") {
        rows.sort((a, b) =>
            String(values[a]).localeCompare(String(values[b]), undefined, { sensitivity: "base", numeric: true }),
        );
    }
    return rows;
}

export type DotPlotLayoutInput = {
    width: number;
    height: number;
    /** Space for the row labels, left of the first panel. */
    yAxisWidth: number;
    /** Space for the field labels, below the panels. */
    xAxisHeight: number;
    /** Space for panel titles, above the panels; 0 when there are none. */
    titleHeight: number;
    panelGap: number;
    panelCount: number;
    fieldCount: number;
    rowCount: number;
};

export type DotPlotLayout = {
    /** Top of the plot area. */
    top: number;
    /** Left edge of each panel's plot area. */
    panelLefts: number[];
    panelWidth: number;
    plotHeight: number;
    cellWidth: number;
    cellHeight: number;
    /** Radius of a dot for the largest fraction. */
    maxRadius: number;
};

export function getDotPlotLayout({
    width,
    height,
    yAxisWidth,
    xAxisHeight,
    titleHeight,
    panelGap,
    panelCount,
    fieldCount,
    rowCount,
}: DotPlotLayoutInput): DotPlotLayout {
    const top = DOT_PLOT_TOP_MARGIN + titleHeight;
    const plotHeight = Math.max(0, height - top - xAxisHeight);
    const gaps = panelGap * Math.max(0, panelCount - 1);
    const panelWidth = panelCount > 0 ? Math.max(0, width - yAxisWidth - DOT_PLOT_RIGHT_MARGIN - gaps) / panelCount : 0;
    const panelLefts = Array.from({ length: panelCount }, (_, p) => yAxisWidth + p * (panelWidth + panelGap));
    const cellWidth = fieldCount > 0 ? panelWidth / fieldCount : 0;
    const cellHeight = rowCount > 0 ? plotHeight / rowCount : 0;
    const maxRadius = Math.min(Math.min(cellWidth, cellHeight) / 2, DOT_PLOT_MAX_RADIUS);
    return { top, panelLefts, panelWidth, plotHeight, cellWidth, cellHeight, maxRadius };
}

/** Dot radius for a fraction, on a square-root scale so that dot area follows the fraction. */
export function getDotRadius(frac: number, fracMax: number, maxRadius: number) {
    if (fracMax <= 0 || frac <= 0) return 0;
    return maxRadius * Math.sqrt(Math.min(frac, fracMax) / fracMax);
}
