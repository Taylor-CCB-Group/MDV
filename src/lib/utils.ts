import type { Param } from "@/charts/ChartTypes";
import type { CategoricalDataType, DataColumn, DataType, GuiSpec, GuiSpecType, GuiValueTypes, LoadedDataColumn, NumberDataType } from "@/charts/charts";
import { type ClassValue, clsx } from "clsx";
import { twMerge } from "tailwind-merge";

export function cn(...inputs: ClassValue[]) {
    return twMerge(clsx(inputs));
}
// nb columnMatchesType has been moved to columnTypeHelpers.ts
/**
 * Check whether given value is an array
 * Acts as a type-predicate, so can be used for narrowing the type of subsequent code.
 */
export function isArray(v: unknown): v is any[] {
    return Array.isArray(v);
}
export function toArray<T>(v: T | T[]) {
    return isArray(v) ? v : [v];
}
//https://github.com/microsoft/TypeScript/issues/45097#issuecomment-882526325
export function notEmpty<TValue>(value: TValue | null | undefined): value is TValue {
    return value !== null && value !== undefined;
}
export type Entries<T> = {
    [K in keyof T]: [K, T[K]];
}[keyof T][];
export function getEntries<T extends object>(o: T) {
    return Object.entries(o) as Entries<T>;
}

export function isDatatypeNumeric(t: DataType): t is NumberDataType {
    return !!t.match(/double|float|int/);
}
export function isDatatypeCategorical(t: DataType): t is CategoricalDataType {
    return t.includes("text");
}
// todo: this is not right...
// export function isGuiValTypeNumeric(t: keyof GuiValueTypes): t is number | [number, number] ? true : false {
//     return !!t.match(/slider|spinner/);
// }

/** concise type-helper factory for making elements of Gui i.e. for settings
 * this will enable type-checking of the GuiSpec objects at the point of declaration.
 * 
 * As of now, this doesn't do anything at runtime - if we want to have some centralised processing of the GuiSpec objects,
 * we could consider adding that here, or may want to keep it separate.
 * 
 * Sometimes the error messages that it generates can be confusing - for example, if you specify
 * a `current_value` for a `"multidropdown"` that is not a `string[]` (which as of this writing is specified
 * as the legal data-type for `GuiSpec<"multidropdown">`), you may get an error message that
 * includes something like
 * ```
 * '"multidropdown"' is assignable to the constraint of type 'T', but 'T' could be instantiated with a different subtype of constraint 'keyof GuiValueTypes'
 * ```
 */
export function g<T extends GuiSpecType>(spec: GuiSpec<T>): GuiSpec<T> {
    return spec;
}
// todo: unit testing
// const m = g({
//     type: "multicolumn",
//     label: "multicolumn test",
//     current_value: ["a"],
//     func(v) {
//         v.map(v => `${v}!`);
//     }
// });

/**
 * Determines whether `str` contains every substring in `substrings`.
 *
 * @param substrings - Array of substrings to look for in `str`
 * @param str - The string to search within
 * @returns `true` if every element of `substrings` is contained in `str`, `false` otherwise
 */
export function stringContainsAll(substrings: string[], str: string) {
    return !substrings.some(i => !str.includes(i));
}

/**
 * Helper function to match two emails, and check if both are same from the start
 */
export function matchEmail(str1: string, str2: string) {
    return str2.includes(str1) && str2.startsWith(str1);
}

/**
 * Parse a delimited string into an array of trimmed, non-empty items.
 *
 * The input is split on commas, spaces, tabs, and newlines; empty segments are removed and each item is trimmed.
 *
 * @param str - The input string containing items separated by commas, spaces, tabs, or newlines
 * @returns An array of trimmed, non-empty strings extracted from `str`
 */
export function parseDelimitedString(str: string) {
    // Split by space or any whitespace character
    const items = str
        .split(/[,\s]+/)
        .map(s => s.trim())
        .filter(Boolean);
    return items;
}

export type NumericRange = [number, number];
export type HistogramScaleMode = "linear" | "log";

type DistributionSummary = {
    lowerTenthShare: number;
    lowerQuarterShare: number;
    lowerHalfShare: number;
    meanPosition: number;
};

/**
 * Decide whether an x-distribution is sufficiently left-heavy to benefit from log display.
 *
 * The thresholds are tuned to prefer `linear` unless the data mass is genuinely concentrated
 * near the low end of the domain.
 */
function shouldUseLogFromSummary(summary: DistributionSummary) {
    return summary.lowerQuarterShare > 0.72 &&
        summary.lowerHalfShare > 0.9 &&
        (summary.lowerTenthShare > 0.35 || summary.meanPosition < 0.22);
}

/**
 * Infer an x-scale mode from raw values by measuring how much probability mass sits in the
 * lower part of the domain.
 *
 * This intentionally avoids choosing `log` based on numeric range alone, so wide-but-even
 * distributions such as coordinates still default to `linear`.
 */
export function resolveAutoHistogramXScaleFromValues(
    domain: NumericRange,
    values: ArrayLike<number> | null | undefined,
): HistogramScaleMode {
    const [min, max] = domain;
    if (!Number.isFinite(min) || !Number.isFinite(max) || min === max) {
        return "linear";
    }
    if (!values || values.length < 32) {
        return "linear";
    }
    let count = 0;
    let lowerTenth = 0;
    let lowerQuarter = 0;
    let lowerHalf = 0;
    let weightedPosition = 0;
    for (let index = 0; index < values.length; index += 1) {
        const rawValue = values[index];
        if (!Number.isFinite(rawValue)) continue;
        const normalized = (rawValue - min) / (max - min);
        if (!Number.isFinite(normalized)) continue;
        const clamped = Math.max(0, Math.min(1, normalized));
        count += 1;
        weightedPosition += clamped;
        if (clamped <= 0.1) lowerTenth += 1;
        if (clamped <= 0.25) lowerQuarter += 1;
        if (clamped <= 0.5) lowerHalf += 1;
    }
    if (count < 32) return "linear";
    return shouldUseLogFromSummary({
        lowerTenthShare: lowerTenth / count,
        lowerQuarterShare: lowerQuarter / count,
        lowerHalfShare: lowerHalf / count,
        meanPosition: weightedPosition / count,
    }) ? "log" : "linear";
}

/**
 * Infer an x-scale mode from histogram bins when only aggregated counts are available.
 *
 * This uses the same left-heavy heuristic as the raw-value version, but estimates the mass
 * position from bin indices instead of original samples.
 */
export function resolveAutoHistogramXScaleFromHistogram(
    domain: NumericRange,
    histogram: number[],
): HistogramScaleMode {
    const [min, max] = domain;
    if (!Number.isFinite(min) || !Number.isFinite(max) || min === max) {
        return "linear";
    }
    const total = histogram.reduce((sum, value) => sum + Math.max(0, value), 0);
    if (total <= 0) return "linear";
    let lowerTenth = 0;
    let lowerQuarter = 0;
    let lowerHalf = 0;
    let weightedPosition = 0;
    histogram.forEach((count, index) => {
        if (count <= 0) return;
        const position = histogram.length <= 1 ? 0.5 : index / (histogram.length - 1);
        weightedPosition += count * position;
        if (position <= 0.1) lowerTenth += count;
        if (position <= 0.25) lowerQuarter += count;
        if (position <= 0.5) lowerHalf += count;
    });
    return shouldUseLogFromSummary({
        lowerTenthShare: lowerTenth / total,
        lowerQuarterShare: lowerQuarter / total,
        lowerHalfShare: lowerHalf / total,
        meanPosition: weightedPosition / total,
    }) ? "log" : "linear";
}

/**
 * Infer whether the y-axis should use log scaling based on how uneven the non-zero bin
 * magnitudes are.
 */
export function resolveAutoHistogramYScale(
    histogram: number[],
): HistogramScaleMode {
    const nonZero = histogram.filter((value) => value > 0);
    if (nonZero.length < 2) return "linear";
    const min = Math.min(...nonZero);
    const max = Math.max(...nonZero);
    const mean = nonZero.reduce((sum, value) => sum + value, 0) / nonZero.length;
    return max / min > 50 || max / Math.max(1, mean) > 10 ? "log" : "linear";
}
