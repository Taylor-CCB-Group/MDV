import type { DataColumn, FieldName } from "@/charts/charts";
import type DataStore from "@/datastore/DataStore";
import type { MulticolumnQuery } from "@/links/link_utils.js";
import type { useState } from "react";
import type { DataType } from "../charts/charts";

// this is more to do with column queries than columns themselves

//new Set(Object.values(BaseChart.types).flatMap(t => t.params).filter(Boolean).flatMap(p => p.type))
// export type Param = "text" | "number" | "multitext" | "text16" | "_multi_column:number" | "_multi_column:all";
export type MultiColumnPrefix = `_multi_column:${"all" | "number"}`; //didn't exist in wild before, but maybe "text" as well?
export type Param = DataType | "number" | MultiColumnPrefix;

// we should be able to describe what permutation of multiple / virtual columns we allow
export type FieldSpec = FieldName | FieldName[] | MulticolumnQuery; //this may be defined elsewhere in future

/** annotation of what kind of column type a given param will accept */
export type CTypes = Param | Param[];
// type MultiPrefix = `_multi_column:${string}`;
type IsMultiParam<T extends CTypes> = T extends Param[] ? true 
: T extends MultiColumnPrefix ? true 
: false;

export type ColumnSelectionProps<T extends CTypes,
    V = IsMultiParam<T> extends true ? FieldSpec : FieldName | MulticolumnQuery> = {
        type: T; //wary of using 'type' as a name - not reserved, but could be confusing
        multiple?: boolean; //also interacts with type "_multi...", perhaps simpler to avoid having both
        setSelectedColumn: (column: V) => void; //what about multiple? also, special values...
        current_value?: V;
        placeholder?: string;
        exclude?: string[];
        dataStore?: DataStore;
    };

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
    props: ColumnSelectionProps<T> & GuiStateProps
): ColumnSelectionProps<T> & GuiStateProps {
    return props;
}
type setBoolean = ReturnType<typeof useState<boolean>>[1];
type GuiStateProps = {
    isExpanded: boolean;
    setIsExpanded: setBoolean;
    isAutocompleteFocused: boolean;
    setIsAutocompleteFocused: setBoolean;
}

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
