// Not necessarily comprehensive set of types, arbitrary coverage of
// - parts of code that intersect with stuff I'm actually using directly in TS
// - general shape of data structures NB - in many cases this can be inferred from JS, so I may remove these again...
//   partly just using them as a form of documentation / notes-to-self.

import type DataStore from "../datastore/DataStore";
import type BaseChart from "./BaseChart";
export type DataType =
    | "integer"
    | "double"
    | "text"
    | "text16"
    | "unique"
    | "multitext"
    | "int32";
export type CategoricalDataType = "text" | "text16" | "multitext";
export type NumberDataType = "integer" | "double" | "int32";

type DataStructureTypes = {
    integer: Uint32Array;
    double: Float32Array; //why is it called 'double'???
    text: Uint8Array;
    text16: Uint16Array;
    multitext: Uint16Array;
    unique: Uint8Array; //raw bytes of strings to be decoded
};
// even if they're just aliases, these could be useful for documentation / clarity
export type ColumnName = string; //this will probably change to ColumnSpecifier with more structured data
// ^^ I'd prefer that to having to reason about the string format everywhere
export type DataSourceName = string;
export type FieldName = string;

type Quantiles = {
    "0.001": [number, number];
    "0.01": [number, number];
    "0.05": [number, number];
};

export type DataColumn<T extends DataType> = {
    name: ColumnName;
    field: FieldName;
    datatype: T;
    data: DataStructureTypes[T];
    values: T extends CategoricalDataType ? string[] : never; //probably wrong for 'unique'
    minMax: T extends NumberDataType ? [number, number] : never;
    quantiles: T extends NumberDataType ? Quantiles : never;
    getValue: (i: number) => T extends CategoricalDataType ? string : number;
};

// export type DataStore = {
//     size: number;
//     filterSize: number;
//     columns: DataColumn<any>[];
//     //tempted to try to be clever about generics in strings (keyof etc)
//     columnGroups: Record<string, {name: string, columns: string[]}>;
//     columnIndex: Record<string, DataColumn<any>>;
//     columnsWithData: string[];
// };

export type RowsAsColumnsLink = {
    name: string;
    name_column: ColumnName;
    subgroups: Record<string, { label: string, name: string, type: "dense" | string }>;
}

export type DataSourceLinks = Record<DataSourceName, {rows_as_columns?: RowsAsColumnsLink}>;

export type DataSource = {
    name: DataSourceName;
    charts: Chart[];
    dataStore: DataStore;
    contentDiv: HTMLDivElement;
    menuBar: HTMLDivElement;
    images?: Record<string, any>;
    regions?: Record<string, any>;
    links?: DataSourceLinks;
};

type DropdownMappedValue<T extends string, V extends string> = {
    [P in T]: string;
} & { [P in V]: string };
/** tuple of `[object[], textKey, valueKey]` where `textKey` and `valueKey` are string properties of the object. */
type DropdownMappedValues<T extends string, V extends string> = [
    Array<DropdownMappedValue<T, V>>,
    T,
    V,
];
/** `'values'` for a dropdown are either a one-element array (with the element being a string array used for both 'text' and 'value'),
 * or a tuple of [object[], textKey, valueKey] where `textKey` and `valueKey` are properties of the object that will be used for 'text'
 * and 'value' respectively.
 */
export type DropDownValues =
    | DropdownMappedValues<string, string>
    | [Array<string>];
export type GuiValueTypes = {
    //is the premise of this correct? technically, could we have a dropdown of numbers?
    dropdown: string;
    multidropdown: string | string[];
    check: boolean;
    text: string;
    textbox: string;
    radiobuttons: string;
    slider: number;
    spinner: number;
    button: undefined;
    doubleslider: [number, number];
    folder: GuiSpec[];
    // color: string; //not a bad idea, copilot... soon...
    // when you are choosing a column, you need to be able to specify
    // - the type of column you can accept
    // - whether you can accept multiple
    // The type you get back ColumnName is going to change into ColumnSpecifier to allow for more complex column references
    // (for now we continue to parse things in a format as in LinkDataDialog)
    // (and we have several keys for different types of column references)
    // There should also be a general way of expressing that a property (like radius) can be set to a 
    // number or a column (with modifiers) - this is where the node editor comes in...
    column: ColumnName;
};
export type GuiSpecType = keyof GuiValueTypes;
export type ColumnSelectionParameters = {
    filter?: DataType[];
    multiple?: boolean;
    exclude?: string[];
}
export type GuiSpec<T extends GuiSpecType> = {
    type: T;
    label: string;
    name: string;
    current_value?: GuiValueTypes[T];
    func?: (v: GuiValueTypes[T]) => void; //optional but ts insists on it?
    values?: T extends "dropdown" | "multidropdown" ? DropDownValues : never;
    // choices is only used for radiobuttons, so we should infer if T is radiobuttons, otherwise never
    choices?: T extends "radiobuttons" ? [string, string][] : never;
    min?: number;
    max?: number;
    step?: number;
    defaultVal?: GuiValueTypes[T];
    columnSelection?: T extends "column" ? ColumnSelectionParameters : never;
};
// todo common interface for AddChartDialog & SettingsDialog - from an end-user perspective, not just types
export type ExtraControl<T extends GuiSpecType> = {
    type: T;
    name: string;
    label: string;
    values?: Array<{ name: string; value: string }>;
    defaultVal?: GuiValueTypes[T];
};
// const a: GuiSpec<'dropdown'> = {
//     type: 'dropdown',
//     label: 'label',
//     name: 'name',
//     current_value: 'current_value',
//     func: (v) => {},
//     values: [[
//         {a: '0a', b: '0b'},
//         {a: '1a', b: '1b'},
//     ], 'a', 'b']
// }
// a.values[0].map(x => x.a);
// ^^ still not well-typed... `x` is `any`.
// biome-ignore lint/suspicious/noRedeclare: I should probably fix this...
interface DataStore {
    getLoadedColumns: () => FieldName[];
    getColumnName: (col: FieldName) => ColumnName | null;
    getColumnList: () => DataColumn[];
    getFilteredIndices: () => Promise<Uint32Array>;
}

export type Chart = {
    getDiv: () => HTMLElement;
    remove: () => void;
    addMenuIcon: (classes: string, info: string) => HTMLElement;
    setSize: (x?: number, y?: number) => void;
    changeBaseDocument: (doc: Document) => void;
    getSettings: () => GuiSpec[];
    removeLayout?: () => void;
    config: any;
    dataStore: DataStore;
    popoutIcon: HTMLElement;
} & BaseChart;

export type ChartState = {
    chart: Chart;
    win?: Window;
    dataSource: DataSource;
};

export type DataSourceSpec = {
    name: DataSourceName;
    dataStore: DataStore;
    //links etc TBD
};

export type ChartManager = {
    charts: Record<string, ChartState>;
    dataSources: DataSourceSpec[];
    dsIndex: Record<string, DataSourceSpec>;
    addMenuIcon: (
        dataSourceName: DataSourceName,
        iconClass: string,
        text: string,
        func: () => void,
    ) => HTMLElement;
    /** probably not something we really want to use publicly like this */
    _popOutChart: (chart: Chart) => void;
};
