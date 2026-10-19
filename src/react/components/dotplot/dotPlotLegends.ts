import { getFractionLegendSize } from "@/react/components/legend/FractionLegend";
import type { FractionLegendSpec } from "@/react/legend/fraction_legend/types";
import { DEFAULT_CONTINUOUS_LEGEND_WIDTH } from "@/react/legend/shared/legendConstants";
import { getContinuousLegendContainerHeight } from "@/react/legend/shared/legendUtils";
import { scaleSqrt } from "d3-scale";
import { getDotRadius } from "./dotPlotLayout";

export const DOT_PLOT_LEGEND_MARGIN = 8;
/** Space `LegendWrapper` adds around a legend's content (padding and border). */
export const LEGEND_WRAPPER_EXTRA = 15;

export type LegendBox = { width: number; height: number };
export type LegendPosition = [left: number, top: number];

/** The fraction legend: dots for a few fractions, the same size as the dots in the panels. */
export function buildFractionLegendSpec(fracMax: number, maxRadius: number): FractionLegendSpec {
    const values =
        fracMax > 0
            ? scaleSqrt()
                  .domain([0, fracMax])
                  .ticks(4)
                  .filter((v) => v > 0)
            : [];
    return {
        label: "% expressing",
        orientation: "horizontal",
        items: values.map((v) => ({ key: String(v), label: String(v), radius: getDotRadius(v, fracMax, maxRadius) })),
    };
}

/** On-screen size of the color legend, as `ColorLegend` sizes a continuous legend. */
export function getColorLegendBox(range: [number, number]): LegendBox {
    return {
        width: DEFAULT_CONTINUOUS_LEGEND_WIDTH + LEGEND_WRAPPER_EXTRA,
        height:
            getContinuousLegendContainerHeight(range, { width: DEFAULT_CONTINUOUS_LEGEND_WIDTH }) +
            LEGEND_WRAPPER_EXTRA,
    };
}

export function getFractionLegendBox(spec: FractionLegendSpec): LegendBox {
    const { width, height } = getFractionLegendSize(spec);
    return { width: width + LEGEND_WRAPPER_EXTRA, height: height + LEGEND_WRAPPER_EXTRA };
}

export type LegendBand = {
    /** Height to keep free above the panels. */
    height: number;
    color: LegendPosition | null;
    fraction: LegendPosition | null;
};

/**
 * Where the legends start: in a band above the panels, the color legend on the left and the
 * fraction legend beside it, or below it when the chart is too narrow for both.
 */
export function getLegendBand(width: number, color: LegendBox | null, fraction: LegendBox | null): LegendBand {
    const m = DOT_PLOT_LEGEND_MARGIN;
    if (color && fraction) {
        if (m + color.width + m + fraction.width + m <= width) {
            return {
                height: m + Math.max(color.height, fraction.height),
                color: [m, m],
                fraction: [m + color.width + m, m],
            };
        }
        return { height: m + color.height + m + fraction.height, color: [m, m], fraction: [m, m + color.height + m] };
    }
    if (color) return { height: m + color.height, color: [m, m], fraction: null };
    if (fraction) return { height: m + fraction.height, color: null, fraction: [m, m] };
    return { height: 0, color: null, fraction: null };
}
