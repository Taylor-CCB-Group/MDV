// Not necessarily comprehensive set of types, arbitrary coverage of
// - parts of code that intersect with stuff I'm actually using directly in TS
// - general shape of data structures NB - in many cases this can be inferred from JS, so I may remove these again...
//   partly just using them as a form of documentation / notes-to-self.

// todo rearrange - maybe have a datastore.d.ts etc
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
    int32: Int32Array;
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
type Colors = (string | number[])[];
type SubgroupName = string;
// wouldn't it be nice to have a type that determins that we know we have data?
export type DataColumn<T extends DataType> = {
    /** human-readable column */
    name: ColumnName; //nb - we should check use of 'name' vs 'field'
    /** id of the column, used internally */
    field: FieldName;
    /** the datatype- can be one of
     *   - `"double"` - any floating point data
     *   - `"integer"` - any integer data
     *   - `"text"` - data containing strings but with no more than 256 categories
     *   - `"text16"` - data containing strings with up to 65536 categories
     *   - `"unique"` - data containing strings but with many categories
     *   - `"multitext"` -
     */
    datatype: T;
    /** whether the column's data can be changed */
    editable?: boolean,
    /**In the case of a double/integer (number) column, the array
     * buffer should be the appropriate size to contain float32s. For text it should be Uint8
     * and contain numbers corresponding to the indexes in the values parameter. For a column of
     * be added later see {@link DataStore#setColumnData}.
     * A {@link LoadedDataColumn<T>} can be used to represent a column that is known to have loaded data.
     */
    data?: DataStructureTypes[T];
    values: T extends CategoricalDataType ? string[] : never; //probably wrong for 'unique'
    /** An array of rgb hex colors. In the case of a `"text"` column the `colors` should match the `values`. For number columns, the list represents
     * colors that will be interpolated. If not specified, default color pallettes will be supplied.
     */
    colors?: Colors;
    /** if `true` then the colors will be displayed on a log scale- useful if the dataset contains outliers.
     * Because a symlog scale is used the data can contain 0 and negative values */
    colorLogScale?: boolean;
    /** the column's values will be displayed as links (text and unique columns only).
     * not sure if this is strictly boolean or can be undefined */
    is_url?: T extends CategoricalDataType ? boolean : never;
    /** the min max values in the column's values (integer/double only) */
    minMax: T extends NumberDataType ? [number, number] : never;
    /** an object describing the 0.05,0.01 and 0,001 qunatile ranges (integer/double only) */
    quantiles: T extends NumberDataType ? Quantiles : never;
    /** if `true` then the store will keep a record that this column has been added and is not permanently stored in the backend */
    dirty?: boolean;
    /** return the value corresponding to a given row index `i`. If the data is categorical, this will be the appropriate value from `values` */
    getValue: (i: number) => T extends CategoricalDataType ? string : number;
    stringLength: T extends "unique" ? number : never;
    delimiter?: T extends "multitext" ? string : never;
    subgroup?: SubgroupName; //not attempting to descriminate the other sg properties being related to this for now
    sgindex?: number;
    sgtype?: "dense" | "sparse"; //?? any other options?
};
export type LoadedDataColumn<T extends DataType> = DataColumn<T> & Required<Pick<DataColumn<T>, "data">>;
function test(column: LoadedDataColumn<"text">) {
    column.data;
    column.values; //why does this degrade to 'any'?
}
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
/** tuple of `[{textKey: string, valueKey: string}[], textKey, valueKey]` where `textKey` and `valueKey` are string properties of the object. */
type DropdownMappedValues<TextKey extends string, ValueKey extends string> = [
    Array<DropdownMappedValue<TextKey, ValueKey>>,
    TextKey,
    ValueKey,
];
/** `'values'` for a dropdown are either a one-element array (with the element being a string array used for both 'text' and 'value'),
 * or a tuple of `[{textKey: string, valueKey: string}[], textKey, valueKey]` where `textKey` and `valueKey` are properties of entries in the array that will be used for 'text'
 * and 'value' respectively.
 */
export type DropDownValues =
    | DropdownMappedValues<string, string>
    | [Array<string>];
export type GuiValueTypes = {
    //is the premise of this correct? technically, could we have a dropdown of numbers?
    dropdown: string;
    multidropdown: string[];
    check: boolean;
    text: string;
    textbox: string;
    radiobuttons: string;
    slider: number;
    spinner: number;
    /** buttons don't have a `current_value` - just the `func()` to call back */
    button: never;
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
    multicolumn: ColumnName[]; //easier to have distinct 'multicolumn' type than overly generic 'column'?
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
    // name: string; //unused - wrongly added in the original type
    current_value: GuiValueTypes[T];
    // is this optional or not? depends slightly how much we lean on mobx current_value mutation for reactivity...
    // should it necessarily not be there for 'folder' - maybe we have some kind of 'is the folder open' property?
    func?: T extends "folder" ? never : (v: GuiValueTypes[T]) => void;
    values?: T extends "dropdown" | "multidropdown" ? DropDownValues : never;
    // choices is only used for radiobuttons, so we should infer if T is radiobuttons, otherwise never
    choices?: T extends "radiobuttons" ? [string, string][] : never;
    //thouht we could make these non-optional and only have it insist for number types, but that's not happening for some reason
    min?: T extends NumberDataType ? number : never;
    max?: T extends NumberDataType ? number : never;
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
/// --- some type-testing code, could be in actual tests...
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
//// biome-ignore lint/suspicious/noRedeclare: I should probably fix this...
// interface DataStore {
//     getLoadedColumns: () => FieldName[];
//     getColumnName: (col: FieldName) => ColumnName | null;
//     getColumnList: () => DataColumn[];
//     getFilteredIndices: () => Promise<Uint32Array>;
//     /** assigns the given `listener` callback to the `id` - will be unceremoniously replaced
//      * if subsequent calls to `addListener` are made with the same `id` */
//     addListener: (id: string, listener: (eventType: string, data: any) => void) => void;
// }

export type Chart<T = any> = {
    getDiv: () => HTMLElement;
    remove: () => void;
    addMenuIcon: (classes: string, info: string) => HTMLElement;
    setSize: (x?: number, y?: number) => void;
    changeBaseDocument: (doc: Document) => void;
    getSettings: () => GuiSpec[];
    removeLayout?: () => void;
    config: T;
    dataStore: DataStore;
    dataSource: DataSource;
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
