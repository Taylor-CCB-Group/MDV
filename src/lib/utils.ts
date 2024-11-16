import type { Param } from "@/charts/ChartTypes";
import type { CategoricalDataType, DataColumn, DataType, GuiValueTypes, NumberDataType } from "@/charts/charts";
import { type ClassValue, clsx } from "clsx";
import { twMerge } from "tailwind-merge";

export function cn(...inputs: ClassValue[]) {
    return twMerge(clsx(inputs));
}
export function columnMatchesType(column: DataColumn<DataType>, type?: Param | Param[]) {
    if (type === undefined) return true;
    if (isArray(type)) return type.some(t => columnMatchesType(column, t));
    if (type === "_multi_column:all") return true;
    // this was allowing a "text" column to be selected for "number" without the '!!'?
    // we should unit-test this function...
    const isNumeric = !!column.datatype.match(/double|float|int/);
    if (type === "_multi_column:number") return isNumeric;
    if (type === "number" && isNumeric) return true;
    return column.datatype === type;
}
/** 
 * Check whether given value {@param v} is an array. 
 * Acts as a type-predicate, so can be used for narrowing the type of subsequent code.
 */
export function isArray(v: unknown): v is any[] {
    return Array.isArray(v);
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
