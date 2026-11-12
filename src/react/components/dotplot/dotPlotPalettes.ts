import { schemeReds } from "d3";

export type DotPlotPalette = "reds" | "blue_white_red";

/** Colour stops for the dot plot's mean expression scale; the same as the legacy `DotPlot.js` uses. */
export const DOT_PLOT_PALETTES: Record<DotPlotPalette, readonly string[]> = {
    reds: schemeReds[8],
    blue_white_red: [
        "#313695",
        "#4575B4",
        "#74ADD1",
        "#ABD9E9",
        "#E0F3F8",
        "#FFFFFF",
        "#FEE090",
        "#FDAE61",
        "#F46D43",
        "#D73027",
        "#A50026",
    ],
};

export function getDotPlotPalette(palette?: string): readonly string[] {
    return palette === "blue_white_red" ? DOT_PLOT_PALETTES.blue_white_red : DOT_PLOT_PALETTES.reds;
}

/** `DataStore.getColumnColors` returns CSS strings, or RGB arrays when asked for them. */
export function toCssColors(colors: readonly (string | number[])[]): string[] {
    return colors.map((c) => (typeof c === "string" ? c : `rgb(${c[0]},${c[1]},${c[2]})`));
}

/** Colour for values that aren't finite, as in `DataStore.getColorFunction`. */
export const MISSING_VALUE_COLOR = "#808080";

/**
 * Maps a value to one of `colors`, binned evenly from `min` to `max`, with the same binning as
 * `DataStore.getColorFunction`. `colors` is the output of `DataStore.getColumnColors` for a
 * numeric column: one colour per bin edge, already interpolated (and log-spaced if asked for).
 */
export function binnedColorScale(colors: readonly (string | number[])[], min: number, max: number) {
    const css = toCssColors(colors);
    const interval = (max - min) / Math.max(1, css.length - 1);
    return (value: number) => {
        if (!Number.isFinite(value)) return MISSING_VALUE_COLOR;
        // e.g. every mean is the same
        if (interval <= Number.EPSILON) return css[0];
        const bin = Math.floor((value - min) / interval);
        return css[Math.min(css.length - 1, Math.max(0, bin))];
    };
}
