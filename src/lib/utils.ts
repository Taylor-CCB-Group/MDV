import type { Param } from "@/charts/ChartTypes";
import type { DataColumn, DataType } from "@/charts/charts";
import { type ClassValue, clsx } from "clsx";
import { twMerge } from "tailwind-merge";

export function cn(...inputs: ClassValue[]) {
    return twMerge(clsx(inputs));
}
export function columnMatchesType(column: DataColumn<DataType>, type?: Param | Param[]) {
    if (type === undefined) return true;
    if (Array.isArray(type)) return type.some(t => columnMatchesType(column, t));
    if (type === "_multi_column:all") return true;
    // this was allowing a "text" column to be selected for "number" without the '!!'?
    // we should unit-test this function...
    const isNumeric = !!column.datatype.match(/double|float|int/);
    if (type === "_multi_column:number") return isNumeric;
    if (type === "number" && isNumeric) return true;
    return column.datatype === type;
}

export function isArray(arr: unknown): arr is any[] {
    return Array.isArray(arr);
}