import type { BaseConfig } from "@/charts/BaseChart";
import type { FieldName, Quantiles } from "@/charts/charts";
import { type FieldSpecs, flattenFields } from "@/lib/columnTypeHelpers";
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
    /** Numeric fields: one cell of the grid per field, or a single plot for one. */
    values: FieldSpecs;
    /**
     * Which prop makes the cells. Several values are always drawn a cell per value (overlaying
     * values × categories in one plot is unreadable), so this only chooses the layout; null is the
     * same as `{ by: "values" }`. `"categories"` (a cell per category, a line per value) is to come.
     */
    grid: MultiLineGridConfig | null;
    /** One x scale for every cell, for values in comparable units (z-scores); otherwise one per cell. */
    share_x: boolean;
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
    /**
     * The brushed range of each value, in data units, by field: each filters the other charts (and
     * the other cells) to the rows with that value in range.
     */
    brush: MultiLineBrushes;
};

export type MultiLineGridConfig = { by: "values" };

/** Brushed ranges by value field; a value with no brush has no entry. */
export type MultiLineBrushes = Partial<Record<FieldName, [number, number]>>;

/** `label` replaces the default axis label (the value's name, or "density"); legacy charts saved one. */
export type MultiLineAxis = AxisConfig & { label: string | null };

/** The legacy chart's axis config, which also saved a `label` (always the column name, or "density"). */
type LegacyAxis = { size?: number; tickfont?: number; textsize?: number; rotate_labels?: boolean; label?: string };

/** Config as saved by the legacy chart or an earlier version of this one. */
export type SavedMultiLineConfig = BaseConfig &
    Partial<Omit<MultiLineConfig, "axis" | "scaletrim" | "categories" | "values" | "brush" | "grid">> & {
        categories?: FieldName | null;
        values?: FieldSpecs;
        /** The legacy chart saved "none" or left it out for no trim. */
        scaletrim?: string | boolean | null;
        axis?: { x?: LegacyAxis; y?: LegacyAxis };
        /** A range for the first value, as the chart saved it when it drew one value; or ranges by field. */
        brush?: [number, number] | Record<string, unknown> | null;
        grid?: { by?: string } | null;
    };

function isRange(value: unknown): value is [number, number] {
    return Array.isArray(value) && value.length === 2 && value.every(Number.isFinite);
}

/** Brushes by field, from a saved config; a single range was the first value's. */
function toBrushes(saved: SavedMultiLineConfig["brush"], values: FieldSpecs): MultiLineBrushes {
    if (isRange(saved)) {
        const field = flattenFields(values)[0];
        return field ? { [field]: saved } : {};
    }
    if (!saved || Array.isArray(saved) || typeof saved !== "object") return {};
    return Object.fromEntries(
        Object.entries(saved).filter((entry): entry is [string, [number, number]] => isRange(entry[1])),
    );
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
    const values = original.values ?? (legacyValue !== undefined ? [legacyValue] : []);
    return {
        ...original,
        param: [],
        categories:
            original.categories !== undefined
                ? original.categories
                : typeof legacyCategory === "string"
                  ? legacyCategory
                  : null,
        values,
        grid: original.grid?.by === "values" ? { by: "values" } : null,
        share_x: original.share_x ?? false,
        stacked,
        fill: original.fill ?? false,
        scaletrim: toScaleTrim(original.scaletrim),
        band_width: typeof original.band_width === "number" && original.band_width > 0 ? original.band_width : null,
        intervals: original.intervals ?? DEFAULT_INTERVALS,
        // the legacy chart hides the legend when stacked, as the rows are labelled
        color_legend: { ...original.color_legend, display: original.color_legend?.display ?? !stacked },
        axis: { x: toAxis(original.axis?.x, 30), y: toAxis(original.axis?.y, 45) },
        fit_labels: original.fit_labels ?? true,
        brush: toBrushes(original.brush, values),
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

export type XDomain = { domain: [number, number]; clamp: boolean };

/**
 * The x domain of each cell, for value ranges in cell order. Shared, every cell gets the union of
 * their domains (the same object), and values are clamped to it when trimmed.
 */
export function getCellXDomains(ranges: readonly ValueRange[], scaletrim: ScaleTrim | null, share: boolean): XDomain[] {
    const own = ranges.map((range) => getXDomain(range, scaletrim));
    if (!share || own.length < 2) return own;
    const shared: XDomain = {
        domain: [Math.min(...own.map((d) => d.domain[0])), Math.max(...own.map((d) => d.domain[1]))],
        clamp: own.some((d) => d.clamp),
    };
    return own.map(() => shared);
}

/** `brushes` with `field`'s brush set to `range`, or removed for null. */
export function withBrush(
    brushes: MultiLineBrushes,
    field: FieldName,
    range: [number, number] | null,
): MultiLineBrushes {
    const { [field]: _, ...rest } = brushes;
    return range ? { ...rest, [field]: range } : rest;
}

/** `brushes` without those of fields not in `fields`. */
export function keepBrushes(brushes: MultiLineBrushes, fields: readonly FieldName[]): MultiLineBrushes {
    return Object.fromEntries(Object.entries(brushes).filter(([field]) => fields.includes(field)));
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

/** Band width slider bounds. It's a fraction of the range whether trimmed or not. */
export const BAND_WIDTH_SLIDER = { min: 0.001, max: 0.5, step: 0.001 } as const;

/** The band width as a fraction of the range, for the slider. */
export function getBandWidthFraction(
    config: Pick<MultiLineConfig, "band_width">,
    domain: [number, number],
    clamp: boolean,
) {
    if (clamp) return config.band_width ?? TRIMMED_DEFAULT_BAND_WIDTH;
    const span = domain[1] - domain[0];
    return span > 0 ? getBandwidth(config, domain, clamp) / span : UNTRIMMED_DEFAULT_BAND_FRACTION;
}

/** `band_width` for a fraction of the range: the fraction itself when trimmed, in data units when not. */
export function bandWidthFromFraction(fraction: number, domain: [number, number], clamp: boolean) {
    const span = domain[1] - domain[0];
    return clamp || !(span > 0) ? fraction : fraction * span;
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
