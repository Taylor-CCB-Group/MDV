import { scaleSqrt } from "d3-scale";
import { getFractionLegendSize } from "@/react/components/legend/FractionLegend";
import type { FractionLegendSpec } from "@/react/legend/fraction_legend/types";
import { DEFAULT_CONTINUOUS_LEGEND_WIDTH } from "@/react/legend/shared/legendConstants";
import { getContinuousLegendContainerHeight } from "@/react/legend/shared/legendUtils";
import { getDotRadius } from "./dotPlotLayout";

export const DOT_PLOT_LEGEND_MARGIN = 8;
/** Space `LegendWrapper` adds around a legend's content (padding and border). */
export const LEGEND_WRAPPER_EXTRA = 15;

export type LegendBox = { width: number; height: number };
export type LegendPosition = [left: number, top: number];

/** How compact the legends are drawn. */
export type LegendSizes = {
    /** Width of the colour legend's content. */
    colorWidth: number;
    /** Most dots the fraction legend shows. */
    fractionItems: number;
};

/** Tried in order until the legends fit in one row: full size first, then more compact. */
const LEGEND_SIZES: LegendSizes[] = [
    { colorWidth: DEFAULT_CONTINUOUS_LEGEND_WIDTH, fractionItems: Number.POSITIVE_INFINITY },
    { colorWidth: DEFAULT_CONTINUOUS_LEGEND_WIDTH, fractionItems: 3 },
    { colorWidth: 150, fractionItems: 3 },
    { colorWidth: 120, fractionItems: 2 },
];

/** Up to `count` of `values`, evenly spread and always keeping the first and last. */
export function pickEvenly<T>(values: T[], count: number): T[] {
    if (values.length <= count) return values;
    if (count <= 1) return values.slice(-1);
    return Array.from({ length: count }, (_, i) => values[Math.round((i * (values.length - 1)) / (count - 1))]);
}

/** The fraction legend: dots for a few fractions, the same size as the dots in the panels. */
export function buildFractionLegendSpec(
    fracMax: number,
    maxRadius: number,
    maxItems = Number.POSITIVE_INFINITY,
): FractionLegendSpec {
    const ticks =
        fracMax > 0
            ? scaleSqrt()
                  .domain([0, fracMax])
                  .ticks(4)
                  .filter((v) => v > 0)
            : [];
    return {
        label: "% expressing",
        orientation: "horizontal",
        items: pickEvenly(ticks, maxItems).map((v) => ({
            key: String(v),
            label: String(v),
            radius: getDotRadius(v, fracMax, maxRadius),
        })),
    };
}

/** On-screen size of the colour legend, as `ColorLegend` sizes a continuous legend. */
export function getColorLegendBox(range: [number, number], width = DEFAULT_CONTINUOUS_LEGEND_WIDTH): LegendBox {
    return {
        width: width + LEGEND_WRAPPER_EXTRA,
        height: getContinuousLegendContainerHeight(range, { width }) + LEGEND_WRAPPER_EXTRA,
    };
}

export function getFractionLegendBox(spec: FractionLegendSpec): LegendBox {
    const { width, height } = getFractionLegendSize(spec);
    return { width: width + LEGEND_WRAPPER_EXTRA, height: height + LEGEND_WRAPPER_EXTRA };
}

export type LegendLayout = LegendSizes & {
    /** Where each legend starts, or null when it isn't shown. */
    color: LegendPosition | null;
    fraction: LegendPosition | null;
};

export type LegendLayoutInput = {
    width: number;
    /** Colour range, or null when the colour legend is hidden. */
    colorRange: [number, number] | null;
    /** Largest fraction, or null when the fraction legend is hidden. */
    fracMax: number | null;
    /** Radius of the largest dot in the panels, so legend dots match them. */
    maxRadius: number;
};

/**
 * Where the legends start over the top left of the plot, as in the legacy dot plot (they can be
 * moved into the chart's strip above the plot instead). Normally in one row, at the largest size
 * that fits the chart's width: the colour legend gets narrower and the fraction legend shows
 * fewer dots before anything is stacked.
 */
export function getLegendLayout({ width, colorRange, fracMax, maxRadius }: LegendLayoutInput): LegendLayout {
    const m = DOT_PLOT_LEGEND_MARGIN;
    const boxesFor = (sizes: LegendSizes) => ({
        color: colorRange ? getColorLegendBox(colorRange, sizes.colorWidth) : null,
        fraction:
            fracMax !== null
                ? getFractionLegendBox(buildFractionLegendSpec(fracMax, maxRadius, sizes.fractionItems))
                : null,
    });
    const fractionLeft = (color: LegendBox | null) => (color ? m + color.width + m : m);

    for (const sizes of LEGEND_SIZES) {
        const { color, fraction } = boxesFor(sizes);
        if (!color && !fraction) return { ...sizes, color: null, fraction: null };
        const rowWidth = fractionLeft(color) + (fraction ? fraction.width + m : 0);
        if (rowWidth > width) continue;
        return {
            ...sizes,
            color: color ? [m, m] : null,
            fraction: fraction ? [fractionLeft(color), m] : null,
        };
    }

    // even the most compact row is too wide: side by side if they fit, else one below the other
    const sizes = LEGEND_SIZES[LEGEND_SIZES.length - 1];
    const { color, fraction } = boxesFor(sizes);
    const sideBySide = !color || !fraction || fractionLeft(color) + fraction.width + m <= width;
    return {
        ...sizes,
        color: color ? [m, m] : null,
        fraction: fraction ? (sideBySide ? [fractionLeft(color), m] : [m, m + (color?.height ?? 0) + m]) : null,
    };
}
