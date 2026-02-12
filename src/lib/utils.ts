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

/**
 * Check whether two strings are exactly equal.
 *
 * @returns `true` if the strings are exactly equal using strict equality (`===`), `false` otherwise.
 */
export function compareStrings(str1: string, str2: string) {
    return str1 === str2;
}