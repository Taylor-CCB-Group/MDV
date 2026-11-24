/**
 * Margins and label text for the multi-line chart's axes. Pure functions; text width is measured
 * through an injected `MeasureText`, so this can be tested without a canvas.
 */
import { scaleLinear } from "d3-scale";
import { type MeasureText, truncateMiddle } from "@/charts/dotplot/xAxisLabels";

/** The font the axes draw with, and that `canvasMeasureText` measures. */
export const AXIS_FONT_FAMILY = "Helvetica, Arial, sans-serif";
export const TICK_LENGTH = 4;
/** Between a tick and its label, and between tick labels and the axis label. */
export const LABEL_GAP = 3;
/** Free space above the plot, and to its right for the last tick label to overhang into. */
export const PAD_TOP = 8;
export const PAD_RIGHT = 12;
/** Text extends this far beyond its font size, for descenders. */
export const LINE_HEIGHT = 1.2;
/** At least this much of the width is left for the plot when fitting row labels. */
const MIN_PLOT_FRACTION = 0.6;
/** Smallest gap between x tick labels. */
const MIN_TICK_SPACING = 10;
/** Most x ticks to try, as for d3's default. */
const MAX_X_TICKS = 10;

/** An axis's margin, used when labels aren't fitted, and its label font size. */
export type AxisSize = { size: number; tickfont: number };

export type MultiLineLayoutInput = {
    width: number;
    height: number;
    axis: { x: AxisSize; y: AxisSize };
    /** Fit margins to the labels; off, the axis `size`s are used. */
    fitLabels: boolean;
    domain: [number, number];
    xLabel: string;
    /** Labels down the left, one per row, top to bottom: when stacked. Otherwise the y label. */
    rowLabels: readonly string[] | null;
    yLabel: string;
};

export type MultiLineLayout = {
    left: number;
    top: number;
    plotWidth: number;
    plotHeight: number;
    /** x tick values, chosen so their labels don't overlap. */
    xTicks: number[];
    xTickFormat: (value: number) => string;
    /** The x axis label, shortened from the middle to the plot's width. */
    xLabel: string;
    /** Row labels shortened to the room on the left, or null when not stacked. */
    rowLabels: string[] | null;
    /** Row label font size: no taller than a row. */
    rowFontSize: number;
    yLabel: string;
};

/** Whether labels centred on `ticks` fit side by side across the plot. */
function ticksFit(
    ticks: number[],
    format: (t: number) => string,
    domain: [number, number],
    plotWidth: number,
    fontSize: number,
    measure: MeasureText,
) {
    if (ticks.length <= 1) return true;
    const widest = Math.max(0, ...ticks.map((t) => measure(format(t), fontSize)));
    // the closest pair of labels, as a fallback may not be evenly spaced
    let gap = Number.POSITIVE_INFINITY;
    for (let i = 1; i < ticks.length; i++) gap = Math.min(gap, ticks[i] - ticks[i - 1]);
    const spacing = (plotWidth * gap) / (domain[1] - domain[0]);
    return widest + MIN_TICK_SPACING <= spacing;
}

/**
 * The most ticks (up to d3's default) whose labels fit side by side across the plot. If even the
 * fewest d3 offers don't fit, just the outer two, then just the one nearest the middle.
 */
function chooseXTicks(domain: [number, number], plotWidth: number, fontSize: number, measure: MeasureText) {
    const scale = scaleLinear().domain(domain);
    for (let count = MAX_X_TICKS; count >= 2; count--) {
        const ticks = scale.ticks(count);
        const format = scale.tickFormat(count);
        if (ticksFit(ticks, format, domain, plotWidth, fontSize, measure)) return { ticks, format };
    }
    const fewest = scale.ticks(2);
    const format = scale.tickFormat(2);
    const outer = fewest.length > 2 ? [fewest[0], fewest[fewest.length - 1]] : fewest;
    if (ticksFit(outer, format, domain, plotWidth, fontSize, measure)) return { ticks: outer, format };
    const middle = (domain[0] + domain[1]) / 2;
    const nearest = fewest.reduce<number[]>(
        (best, t) => (best.length && Math.abs(best[0] - middle) <= Math.abs(t - middle) ? best : [t]),
        [],
    );
    return { ticks: nearest, format };
}

export function getMultiLineLayout(input: MultiLineLayoutInput, measure: MeasureText): MultiLineLayout {
    const { width, height, axis, fitLabels, domain, rowLabels } = input;
    const xFont = axis.x.tickfont;
    const yFont = axis.y.tickfont;

    // bottom: tick, tick labels, gap, axis label
    const bottom = fitLabels ? TICK_LENGTH + LABEL_GAP + 2 * xFont * LINE_HEIGHT + LABEL_GAP : axis.x.size;
    const plotHeight = Math.max(0, height - PAD_TOP - bottom);
    const rowHeight = rowLabels?.length ? plotHeight / rowLabels.length : plotHeight;
    const rowFontSize = Math.max(1, Math.min(yFont, rowHeight));

    // left: row labels when stacked, else the rotated y label
    let left: number;
    let fittedRowLabels: string[] | null = null;
    if (rowLabels) {
        const room = fitLabels
            ? width * (1 - MIN_PLOT_FRACTION) - TICK_LENGTH - LABEL_GAP
            : axis.y.size - TICK_LENGTH - LABEL_GAP;
        fittedRowLabels = rowLabels.map((label) => truncateMiddle(label, Math.max(0, room), rowFontSize, measure));
        const widest = Math.max(0, ...fittedRowLabels.map((label) => measure(label, rowFontSize)));
        left = fitLabels ? Math.ceil(widest + TICK_LENGTH + LABEL_GAP + 1) : axis.y.size;
    } else {
        left = fitLabels ? Math.ceil(yFont * LINE_HEIGHT + LABEL_GAP + 1) : axis.y.size;
    }

    const plotWidth = Math.max(0, width - left - PAD_RIGHT);
    const { ticks, format } = chooseXTicks(domain, plotWidth, xFont, measure);
    return {
        left,
        top: PAD_TOP,
        plotWidth,
        plotHeight,
        xTicks: ticks,
        xTickFormat: format,
        xLabel: truncateMiddle(input.xLabel, plotWidth, xFont, measure),
        rowLabels: fittedRowLabels,
        rowFontSize,
        yLabel: truncateMiddle(input.yLabel, plotHeight, yFont, measure),
    };
}
