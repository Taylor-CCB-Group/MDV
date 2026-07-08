/**
 * Color-legend-specific spec types used by the generic React legend wrapper.
 *
 * These types are intentionally scoped to the color legend path for phase 1.
 * Future legend migrations (fraction/node/link/etc.) can add sibling spec types
 * and reuse the same wrapper pattern without changing existing color behavior.
 */
export type ColorLegendCategoricalItem = {
    color: string;
    name: string;
    value: string;
};

export type ColorLegendCategoricalSpec = {
    kind: "categorical";
    label: string;
    column: string;
    items: ColorLegendCategoricalItem[];
};

export type ColorLegendContinuousSpec = {
    kind: "continuous";
    label: string;
    column: string;
    colors: string[];
    range: [number, number];
    width?: number;
    height?: number;
};

export type ColorLegendSpec =
    | ColorLegendCategoricalSpec
    | ColorLegendContinuousSpec;

export type ColorLegendOverrideValues = {
    min?: number | null;
    max?: number | null;
    colorLogScale?: boolean;
    colors?: string[];
    bins?: number;
};

export type ColorLegendBuildConfig = {
    name?: string;
    overrideValues?: ColorLegendOverrideValues;
    /** @deprecated Use overrideValues instead. Kept for existing color config callers. */
    overideValues?: ColorLegendOverrideValues;
};
