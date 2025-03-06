import type { DataColumn, FieldName, LoadedDataColumn } from "@/charts/charts";
import type DataStore from "@/datastore/DataStore";
import type { MultiColumnQuery } from "@/links/link_utils.js";
import type { useState } from "react";
import type { DataType } from "../charts/charts";
import { isArray } from "./utils";

// this is more to do with column queries than columns themselves

//new Set(Object.values(BaseChart.types).flatMap(t => t.params).filter(Boolean).flatMap(p => p.type))
// export type Param = "text" | "number" | "multitext" | "text16" | "_multi_column:number" | "_multi_column:all";
export type MultiColumnPrefix = `_multi_column:${"all" | "number"}`; //didn't exist in wild before, but maybe "text" as well?
export type Param = DataType | "number" | MultiColumnPrefix;

// we should be able to describe what permutation of multiple / virtual columns we allow
// the assumption that if we have an array it will only be strings is wrong.
// up for review & documentation...
export type FieldSpecs = (MultiColumnQuery | FieldName)[];
//@ts-expect- error we need to be more basic & explicit about things that are and are not arrays (should have SingleColumnQuery)
export type FieldSpec = FieldName | MultiColumnQuery;// | (FieldName | MultiColumnQuery)[];
/**
 * Should the type of spec be more restricted...
 * @returns an array of strings representing 'field names' - column IDs.
 */
export function flattenFields(spec: FieldSpec | FieldSpecs): FieldName[] {
    return isArray(spec) ? spec.flatMap(flattenFields) : typeof spec === "string" ? [spec] : spec.fields;
}
/** annotation of what kind of column type a given param will accept */
export type CTypes = Param | Param[];
// type MultiPrefix = `_multi_column:${string}`;
export type IsMultiParam<T extends CTypes> = T extends Param[] ? true 
: T extends MultiColumnPrefix ? true 
: false;

export type ColumnSelectionProps<T extends CTypes,
    V = IsMultiParam<T> extends true ? FieldSpecs : FieldName | MultiColumnQuery> = { //should we have a SingleColumnQuery?
        type?: T; //wary of using 'type' as a name - not reserved, but could be confusing. also wary of optional type
        multiple?: boolean; //also interacts with type "_multi...", perhaps simpler to avoid having both
        setSelectedColumn: (column: V) => void;
        current_value?: V;
        placeholder?: string;
        exclude?: string[];
        dataStore?: DataStore;
    };

/** 
 * Given a `ColumnSelectionProps` object, this function should return the same object with the correct type annotations.
 */
export function inferGenericColumnSelectionProps<T extends CTypes>(
    props: ColumnSelectionProps<T>
): ColumnSelectionProps<T> {
    return props;
}

/// these result in the setSelectedColumn function being typed correctly
// const testMulti = createGenericColumnSelectionProps({
//     type: "_multi_column:number",
//     setSelectedColumn: (multiColumns) => { },
// });
// const testSingle = createGenericColumnSelectionProps({
//     type: "number",
//     setSelectedColumn: (column) => { },
// });
// const testArray = createGenericColumnSelectionProps({
//     type: ["number", "text"],
//     setSelectedColumn: (multiColumns) => { },
// });

export function isMultiColumn(type: CTypes): type is MultiColumnPrefix | Param[] {
    return typeof type === "string" && type.startsWith("_multi_column:");
}
export function inferGenericColumnGuiProps<T extends CTypes>(
    props: ColumnSelectionProps<T>
): ColumnSelectionProps<T> {
    return props;
}
/**
 * This is for filtering columns based on some relatively complex type specification potentially including things like `"_multi_column:number"`...
 * as of this writing, it is not able to function as a type-predicate. For more type-narrowing in more simple cases where you want to assert
 * that a column is of a particular @link {DataType}, @see {isColumnOfType}
 */
export function columnMatchesType(column: DataColumn<DataType>, type?: Param | Param[]): boolean {
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
/**
 * Checks whether @param column has a `datatype` precisely matching the given @param datatype
 */
export function isColumnOfType<T extends DataType>(column: DataColumn<DataType>, datatype: T): column is DataColumn<T> {
    return column.datatype === datatype;
}
/**
 * Checks whether the given @param column exists and has `data`, also acting as a type-predicate such that
 * subsequent references to that column can safely use it.
 */
export function isColumnLoaded(column?: DataColumn<DataType>): column is LoadedDataColumn<DataType> {
    return column?.data !== undefined;
}
export function allColumnsLoaded(columns: DataColumn<DataType>[]): columns is LoadedDataColumn<DataType>[] {
    return columns.every(isColumnLoaded);
}

