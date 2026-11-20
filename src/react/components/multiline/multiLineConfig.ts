import type { BaseConfig } from "@/charts/BaseChart";
import type { FieldName, Quantiles } from "@/charts/charts";
import type { FieldSpecs } from "@/lib/columnTypeHelpers";
import type { AxisConfig } from "@/react/scatter_state";

/** The legacy `multi_line_chart` converts to this (see `legacyCharts.ts`). */
export const MULTI_LINE_REACT_TYPE = "multi_line_chart_react";

/** A percentile to trim the x axis to, as a key of a column's `quantiles`. */
export type ScaleTrim = keyof Quantiles;
export const SCALE_TRIMS: readonly ScaleTrim[] = ["0.001", "0.01", "0.05"];

/** Bandwidth when trimmed, as a fraction of the trimmed range: the legacy chart's default. */
const TRIMMED_DEFAULT_BAND_WIDTH = 0.1;
/** Bandwidth when not trimmed, as a fraction of the column's range: the legacy chart's default. */
const UNTRIMMED_DEFAULT_BAND_FRACTION = 0.01;
const DEFAULT_INTERVALS = 40;

/**
 * Density curves of numeric values, one per category. Uses named props rather than `param`,
 * which stays empty.
 */
export type MultiLineConfig = BaseConfig & {
    /** Category column: one line per value of it. */
    categories: FieldName | null;
    /** Numeric fields. One is drawn for now; several, as a grid, come later. */
    values: FieldSpecs;
    /** One row per category, each scaled to its own height, instead of overlaid lines. */
    stacked: boolean;
    /** Fill under each line. */
    fill: boolean;
    /** Trim the x axis to this percentile of the values; null for the column's full range. */
    scaletrim: ScaleTrim | null;
    /**
     * Half-width of the kernel: in data units, or as a fraction of the trimmed range when
     * `scaletrim` is set (as in the legacy chart). Null for the default.
     */
    band_width: number | null;
    /** Roughly how many points each line has (a hint for d3's `ticks`). */
    intervals: number;
    color_legend: { display: boolean; pos?: [number, number] };
    axis: { x: MultiLineAxis; y: MultiLineAxis };
    /** Fit the axis margins to the labels; off, the axis sizes are used, as in the legacy chart. */
    fit_labels: boolean;
    /** The brushed range of the value, in data units, which filters other charts; null for none. */
    brush: [number, number] | null;
};

/** `label` replaces the default axis label (the value's name, or "density"); legacy charts saved one. */
export type MultiLineAxis = AxisConfig & { label: string | null };

/** The legacy chart's axis config, which also saved a `label` (always the column name, or "density"). */
type LegacyAxis = { size?: number; tickfont?: number; textsize?: number; rotate_labels?: boolean; label?: string };

/** Config as saved by the legacy chart or an earlier version of this one. */
export type SavedMultiLineConfig = BaseConfig &
    Partial<Omit<MultiLineConfig, "axis" | "scaletrim" | "categories" | "values">> & {
        categories?: FieldName | null;
        values?: FieldSpecs;
        /** The legacy chart saved "none" or left it out for no trim. */
        scaletrim?: string | boolean | null;
        axis?: { x?: LegacyAxis; y?: LegacyAxis };
    };

function isRange(value: unknown): value is [number, number] {
    return Array.isArray(value) && value.length === 2 && value.every(Number.isFinite);
}

function toScaleTrim(value: unknown): ScaleTrim | null {
    return SCALE_TRIMS.find((t) => t === value) ?? null;
}

function toAxis(axis: LegacyAxis | undefined, size: number) {
    return {
        size: axis?.size ?? size,
        tickfont: axis?.tickfont ?? axis?.textsize ?? 13,
        rotate_labels: axis?.rotate_labels ?? false,
        label: axis?.label || null,
    };
}

/**
 * Named props from a saved config. A legacy `param` is `[value, category]`, which is also what
 * the add-chart dialog builds from this chart type's `params`.
 */
export function adaptMultiLineConfig(original: SavedMultiLineConfig): MultiLineConfig {
    const param = Array.isArray(original.param) ? original.param : [];
    const legacyValue = param[0];
    const legacyCategory = param[1];
    const stacked = original.stacked ?? false;
    return {
        ...original,
        param: [],
        categories:
            original.categories !== undefined
                ? original.categories
                : typeof legacyCategory === "string"
                  ? legacyCategory
                  : null,
        values: original.values ?? (legacyValue !== undefined ? [legacyValue] : []),
        stacked,
        fill: original.fill ?? false,
        scaletrim: toScaleTrim(original.scaletrim),
        band_width: typeof original.band_width === "number" && original.band_width > 0 ? original.band_width : null,
        intervals: original.intervals ?? DEFAULT_INTERVALS,
        // the legacy chart hides the legend when stacked, as the rows are labelled
        color_legend: { ...original.color_legend, display: original.color_legend?.display ?? !stacked },
        axis: { x: toAxis(original.axis?.x, 30), y: toAxis(original.axis?.y, 45) },
        fit_labels: original.fit_labels ?? true,
        brush: isRange(original.brush) ? original.brush : null,
    };
}

/** What a value column offers for the x axis, from its metadata: no data needs to be loaded. */
export type ValueRange = {
    minMax: [number, number];
    quantiles?: Partial<Record<ScaleTrim, [number, number]>>;
};

/** The x domain to draw over, in data units, and whether the values are clamped to it. */
export function getXDomain(range: ValueRange, scaletrim: ScaleTrim | null) {
    const trimmed = scaletrim ? range.quantiles?.[scaletrim] : undefined;
    if (trimmed && Number.isFinite(trimmed[0]) && Number.isFinite(trimmed[1]) && trimmed[1] > trimmed[0]) {
        return { domain: trimmed, clamp: true };
    }
    return { domain: range.minMax, clamp: false };
}

/**
 * Kernel half-width in data units. `band_width` is a fraction of the trimmed range when trimmed
 * (so the same setting suits any column), and in data units otherwise, as in the legacy chart.
 */
export function getBandwidth(config: Pick<MultiLineConfig, "band_width">, domain: [number, number], clamp: boolean) {
    const span = domain[1] - domain[0];
    if (clamp) return (config.band_width ?? TRIMMED_DEFAULT_BAND_WIDTH) * span;
    return config.band_width ?? span * UNTRIMMED_DEFAULT_BAND_FRACTION;
}

/** Band width slider bounds and the value it shows, in the units `band_width` is in. */
export function getBandWidthSetting(
    config: Pick<MultiLineConfig, "band_width">,
    domain: [number, number],
    clamp: boolean,
) {
    if (clamp) {
        return { min: 0.005, max: 0.5, step: 0.005, value: config.band_width ?? TRIMMED_DEFAULT_BAND_WIDTH };
    }
    const span = domain[1] - domain[0] || 1;
    return { min: span / 1000, max: span / 10, step: span / 1000, value: getBandwidth(config, domain, clamp) };
}

/**
 * The rows a brush selects, as a value range. Values beyond a trimmed axis are drawn at its ends,
 * so a brush that reaches an end takes in everything beyond it; so does one reaching the end of an
 * untrimmed axis, where the brush can stop a rounding error short of the last value.
 */
export function getBrushFilterRange(brush: [number, number], domain: [number, number]): { min: number; max: number } {
    const lo = Math.min(brush[0], brush[1]);
    const hi = Math.max(brush[0], brush[1]);
    const tolerance = (domain[1] - domain[0]) * 1e-6;
    return {
        min: lo <= domain[0] + tolerance ? Number.NEGATIVE_INFINITY : lo,
        max: hi >= domain[1] - tolerance ? Number.POSITIVE_INFINITY : hi,
    };
}
