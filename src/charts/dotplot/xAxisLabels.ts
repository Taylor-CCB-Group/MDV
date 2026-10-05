/**
 * X-axis label layout for the dot plot: which text to show, at what angle and
 * font size, and how much bottom margin it needs. Pure functions; text width
 * is measured through an injected `MeasureText` so this can be tested without
 * a canvas.
 */

export type MeasureText = (text: string, fontSize: number) => number;

export type XLabelLayoutInput = {
    /** labels in column order */
    labels: string[];
    /** left margin of the plot area, i.e. how far the first column is from the chart edge */
    plotLeft: number;
    plotWidth: number;
    chartHeight: number;
    /** the user's x axis text size; labels never grow past it */
    fontSize: number;
    /** space for the x axis title, 0 when there is none */
    titleSize: number;
    /** shorten long labels and cap the margin */
    truncate: boolean;
    /** label only every nth column when they would overlap */
    thin: boolean;
};

export type XLabelLayout = {
    angle: 0 | 45;
    fontSize: number;
    /** bottom margin needed for the labels */
    margin: number;
    /** labels to display, in column order; "" for unlabelled columns */
    labels: string[];
};

// d3 axisBottom draws tick text this far below the axis line (tick size 6 + padding 3)
export const X_TICK_OFFSET = 9;
export const X_LABEL_GAP = 4;
// bottom margin for x labels is capped at the smaller of these
export const X_MARGIN_MAX_PX = 90;
export const X_MARGIN_MAX_FRACTION = 0.3;
// with truncation off, the margin grows to fit full labels up to this fraction
export const X_MARGIN_UNTRUNCATED_MAX_FRACTION = 0.6;
// x labels shrink to fit narrow columns, down to this size
export const X_MIN_FONT_SIZE = 7;
// gap between neighbouring 45° labels, as a multiple of font size
export const X_LABEL_SPACING = 1.1;

let measureContext: CanvasRenderingContext2D | null = null;

/** Measures text as the axis draws it (Helvetica at the given size). */
export const canvasMeasureText: MeasureText = (text, fontSize) => {
    measureContext ??= document.createElement("canvas").getContext("2d");
    if (!measureContext) {
        return text.length * fontSize * 0.6;
    }
    measureContext.font = `${fontSize}px Helvetica`;
    return measureContext.measureText(text).width;
};

export type XAxisColumn = {
    field: string;
    name: string;
    /** whether the column belongs to a rows-as-columns subgroup */
    inSubgroup: boolean;
};

/**
 * Column names for the x axis. When every column comes from the same
 * subgroup, the repeated `(subgroup)` suffix is dropped to save space.
 */
export function stripSharedSubgroupSuffix(columns: XAxisColumn[]): string[] {
    const names = columns.map((c) => c.name);
    const keys = new Set(
        columns.map((c) => (c.inSubgroup ? c.field.split("|")[0] : null)),
    );
    if (keys.size !== 1) {
        return names;
    }
    const [key] = keys;
    if (key === null || key === undefined) {
        return names;
    }
    const suffix = `(${key})`;
    return names.map((name) =>
        name.endsWith(suffix) && name.length > suffix.length
            ? name.slice(0, -suffix.length).trimEnd()
            : name,
    );
}

/**
 * Shortens text to fit maxWidth by cutting from the middle, so labels that
 * share a prefix (e.g. Ensembl ids) stay distinguishable. Binary-searches the
 * number of characters kept, since width grows with each character added.
 */
export function truncateMiddle(
    text: string,
    maxWidth: number,
    fontSize: number,
    measure: MeasureText,
): string {
    if (measure(text, fontSize) <= maxWidth) {
        return text;
    }
    const truncate = (keep: number) =>
        `${text.slice(0, Math.ceil(keep / 2))}…${text.slice(text.length - Math.floor(keep / 2))}`;
    let lo = 0;
    let hi = text.length - 1;
    while (lo < hi) {
        const mid = Math.ceil((lo + hi) / 2);
        if (measure(truncate(mid), fontSize) <= maxWidth) {
            lo = mid;
        } else {
            hi = mid - 1;
        }
    }
    return truncate(lo);
}

/**
 * Lays x labels flat when they all fit within a column, otherwise at 45°.
 * At 45° the font shrinks to fit the column width (down to a minimum); with
 * `thin`, only every nth column is labelled below that. With `truncate`, the
 * margin is capped and labels past the cap, or reaching past the left edge of
 * the chart, are shortened from the middle; otherwise the margin grows to fit
 * the full labels.
 */
export function layoutXLabels(
    input: XLabelLayoutInput,
    measure: MeasureText,
): XLabelLayout {
    const { labels, plotLeft, chartHeight, titleSize, truncate } = input;
    const bandWidth = input.plotWidth / labels.length;

    const baseWidths = labels.map((l) => measure(l, input.fontSize));
    if (Math.max(...baseWidths) <= bandWidth - X_LABEL_GAP) {
        return {
            angle: 0,
            fontSize: input.fontSize,
            margin: Math.ceil(X_TICK_OFFSET + input.fontSize + titleSize + X_LABEL_GAP),
            labels,
        };
    }

    // perpendicular distance between neighbouring 45° labels
    const spacing = bandWidth * Math.SQRT1_2;
    const fontSize = Math.max(
        X_MIN_FONT_SIZE,
        Math.min(input.fontSize, spacing / X_LABEL_SPACING),
    );
    const step = input.thin
        ? Math.ceil((fontSize * X_LABEL_SPACING) / spacing)
        : 1;

    const fixed = X_TICK_OFFSET + fontSize * Math.SQRT1_2 + titleSize + X_LABEL_GAP;
    const cap = Math.max(
        fixed + fontSize,
        Math.min(X_MARGIN_MAX_PX, Math.floor(chartHeight * X_MARGIN_MAX_FRACTION)),
    );
    const capWidth = (cap - fixed) / Math.SQRT1_2;
    let longest = 0;
    const display = labels.map((label, i) => {
        if (i % step !== 0) {
            return "";
        }
        // at 45° a label reaches left from its tick by width·cos45
        const tickX = plotLeft + (i + 0.5) * bandWidth;
        const maxTextWidth = truncate
            ? Math.min(capWidth, (tickX - X_LABEL_GAP) / Math.SQRT1_2)
            : Number.POSITIVE_INFINITY;
        const text = truncateMiddle(label, maxTextWidth, fontSize, measure);
        longest = Math.max(longest, measure(text, fontSize));
        return text;
    });
    // without truncation, still leave room for the dots
    const margin = Math.min(
        fixed + longest * Math.SQRT1_2,
        chartHeight * X_MARGIN_UNTRUNCATED_MAX_FRACTION,
    );
    return { angle: 45, fontSize, margin: Math.ceil(margin), labels: display };
}
