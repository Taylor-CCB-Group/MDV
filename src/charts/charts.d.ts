// Not necessarily comprehensive set of types, arbitrary coverage of
// - parts of code that intersect with stuff I'm actually using directly in TS
// - general shape of data structures NB - in many cases this can be inferred from JS, so I may remove these again...
//   partly just using them as a form of documentation / notes-to-self.

// todo rearrange - maybe have a datastore.d.ts etc
import type { CTypes, FieldSpec, FieldSpecs } from "@/lib/columnTypeHelpers";
import type DataStore from "../datastore/DataStore";
import type BaseChart from "./BaseChart";
// import type { DataType } from "../datatypes";
/**
 * The are the names used to refer to the types of data can be stored in a column.
 */
export type DataType = 
    | "integer"
    | "double"
    | "text"
    | "text16"
    | "unique"
    | "multitext"
    | "int32";
/**
 * To denote cases where any categorical data can be used,
 * this is the union of the {@link DataType}s that share related properties.
 */
export type CategoricalDataType = "text" | "text16" | "multitext";
/**
 * To denote cases where any numerical data can be used,
 * this is the union of the {@link DataType}s that share related properties.
 */
export type NumberDataType = "integer" | "double" | "int32";

/**
 * Associates a {@link DataType} with the corresponding `TypedArray` type that will be used to store the data.
 * 
 * This should be associated with definitions in `datatypes.ts`.
 */
type DataStructureTypes = {
    int32: Int32Array;
    double: Float32Array; //why is it called 'double'???
    integer: Uint32Array;
    text16: Uint16Array;
    text: Uint8Array;
    unique: Uint8Array; //raw bytes of strings to be decoded
    multitext: Uint8Array;
};
export type { DataType, DataStructureTypes };
/**
 * Union of all the possible data structures that can be used to store column data.
 */
export type ColumnData = DataStructureTypes[DataType];
// even if they're just aliases, these could be useful for documentation / clarity
export type ColumnName = string; //this will probably change to ColumnSpecifier with more structured data
// ^^ I'd prefer that to having to reason about the string format everywhere
export type DataSourceName = string;
export type FieldName = string;

export type Quantiles = {
    "0.001": [number, number];
    "0.01": [number, number];
    "0.05": [number, number];
};
type Colors = (string | number[])[];
type SubgroupName = string;
/**
 * Represents a column of data in a {@link DataStore}.
 */
export type DataColumn<T extends DataType> = {
    /** human-readable column name, to be displayed in GUI etc */
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
    editable?: boolean;
    /**In the case of a double/integer (number) column, the array
     * buffer should be the appropriate size to contain float32s. For text it should be Uint8
     * and contain numbers corresponding to the indexes in the values parameter. For a column of
     * type unique it should be a JavaScript array. This parameter is optional as the data can
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
export type LoadedDataColumn<T extends DataType> = DataColumn<T> &
    Required<Pick<DataColumn<T>, "data">>;
// function test(column: LoadedDataColumn<"text">) {
//     column.data;
//     column.values; //why does this degrade to 'any'?
// }
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
    name_column: FieldName;
    subgroups: Record<
        string,
        { label: string; name: string; type: "dense" | string }
    >;
};

export type DataSourceLinks = Record<
    DataSourceName,
    { rows_as_columns?: RowsAsColumnsLink }
>;

/**
 * This is a very provisional type for a way of specifying a zarr store that can be used to load data.
 * 
 * I would prefer the definition for this to be in a schema, will need to get back to this later.
 */
export type ExperimentalZarrStore = {
    path: string,
    type: "xe_cells" | "xe_transcripts" | "xe_analysis"
}

/**
 * nb - there should be a type for the entries in `datasources.json`...
 * *but this is not that type!* even though it is being used in places where that type is expected.
 * There is some confusion here about the internal representation in `ChartManager` of a panel with
 * charts, menuBar etc, and the representation of a datasource in `datasources.json`.
 * We may add a zod schema for the latter, and re-evaluate where this is used.
 */
export type DataSource = {
    name: DataSourceName;
    charts: Chart[]; //what's this? probably an artifact of me not understanding when I first wrote this...
    dataStore: DataStore;
    contentDiv: HTMLDivElement;
    menuBar: HTMLDivElement;
    images?: Record<string, any>;
    regions?: Record<string, any>;
    links?: DataSourceLinks;
    size: number;
    columns: DataColumn<DataType>[];
};// | ExperimentalZarrStore; ? something something spatialdata.js ...
// maybe that the

type DropdownMappedValue<T extends string, V extends string> = {
    [P in T]: string;
} & { [P in V]: string };
/** tuple of `[{textKey: string, valueKey: string}[], textKey, valueKey]` where `textKey` and `valueKey` are string properties of the object. */
export type DropdownMappedValues<TextKey extends string, ValueKey extends string> = [
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
    // The type you get back FieldName is going to change into ColumnSpecifier to allow for more complex column references
    // (for now we continue to parse things in a format as in LinkDataDialog)
    // (and we have several keys for different types of column references)
    // There should also be a general way of expressing that a property (like radius) can be set to a
    // number or a column (with modifiers) - this is where the node editor comes in...
    column: FieldSpec;
    multicolumn: FieldSpecs; //easier to have distinct 'multicolumn' type than overly generic 'column'?
};
export type GuiSpecType = keyof GuiValueTypes;
type GV<T extends GuiSpecType> = GuiValueTypes[T];
type GuiFunc<T extends GuiSpecType> = (v: GV<T>) => void | Promise<void>;
// type GuiSpecExperiment<T extends keyof GuiSpecType> = T extends infer K
//     ? K extends keyof GuiValueTypes
//         ? {
//             readonly type: K;
//             label: string;
//             current_value: GuiValueTypes[K];
//             func?: GuiFunc<K>;
//         }
//         : never
//     : never;

/**
 * This type describes an object that can be used to specify a GUI element in the settings dialog.
 * The type parameter `T` is used to specify the type of the GUI element (such as `"slider"`, `"dropdown"` etc),
 * from which we can infer `GV<T>` (`GuiValueType[T]`) used to specify the type of the `current_value` property of the object,
 * as well as some other associated properties that might be expected.
 * 
 * For example, a `GuiSpec<"slider">` has `current_value: number`,
 * in which case (because `GV<T> extends number`) we may also expect to have associated properties like `min` & `max`.
 * The `GuiValueType["doubleslider"]` is `[number, number]`, which also means that there are `min` & `max` properties etc.
 * For other types, `min` & `max` are not expected, so they are `never`.
 * 
 * There is a type-helper function `g({...})` that can be used when creating these objects, which can help with type inference,
 * and allow for editor feedback to suggest properties etc.
 */
export type GuiSpec<T extends GuiSpecType> = {
    readonly type: T;
    label: string;
    current_value: GV<T>;
    // is this optional or not? depends slightly how much we lean on mobx current_value mutation for reactivity...
    func?: GuiFunc<T>;
    //@ts-check !this is for review... it should *not* be optional, but if I make it non-optional then it demands values for 'never'...
    values?: T extends "dropdown" | "multidropdown" ? DropDownValues : never;
    /**
     * Used by `"radiobuttons"` to specify the choices available.
     * This is a tuple of `[string, string][]` where the first element of each tuple is the label to display,
     * and the second element is the value to be used when the radio button is selected.
     */
    choices?: T extends "radiobuttons" ? [string, string][] : never;
    //thouht we could make these non-optional and only have it insist for number types, but that's not happening for some reason
    min?: GV<T> extends number
        ? number
        : GV<T> extends [number, number]
          ? number
          : never;
    max?: GV<T> extends number
        ? number
        : GV<T> extends [number, number]
          ? number
          : never;
    step?: GV<T> extends number
        ? number
        : GV<T> extends [number, number]
          ? number
          : never;
    continuous?: GV<T> extends number
        ? boolean
        : GV<T> extends [number, number]
          ? boolean
          : never;
    defaultVal?: GV<T>;
    /**
     * Expresses the type of column that can be selected in a column selection dialog.
     * Can be in the form of {@link DataType} as in {@link DataColumn#datatype}, 
     * or other form like `"_multi_column:number" | "number"`.
     */
    columnType?: T extends ("column" | "multicolumn") ? CTypes : never;
};
// export type GuiSpecs = Array<GuiSpec<GuiSpecType>>;
/**
 * This creates a discriminated union of `GuiSpec<"folder"> | GuiSpec<"slider"> | ...`, which can be preferable to the more
 * obvious seeming `GuiSpec<GuiSpecType>`, which evaluates to `GuiSpec<"folder" | "slider" | ...>` - causing
 * some obscure co/contra-variance issues, often resulting in irritating and incomprehensible error messages.
 * 
 * For example, in the Settings Dialog, we have an array of settings - which could be any kind of gui spec - and
 * we want to iterate over them and pass them to `AbstractComponent` - which will then render the appropriate
 * component for each setting, while also inferring that the types are correct.
 * 
 * If we have `GuiSpec<GuiSpecType>` then we can't pass that to `AbstractComponent` because it will complain,
 * with a lengthy error message that includes something like
 * ```
 * Type 'AnyGuiSpec' is not assignable to type 'GuiSpec<keyof GuiValueTypes>'.
 *  Type 'GuiSpec<"text">' is not assignable to type 'GuiSpec<keyof GuiValueTypes>'.
 *      Types of property 'func' are incompatible.
 *      Type 'GuiFunc<"text"> | undefined' is not assignable to type 'GuiFunc<keyof GuiValueTypes> | undefined'.
 *          Type 'GuiFunc<"text">' is not assignable to type 'GuiFunc<keyof GuiValueTypes>'.
 *          Type 'keyof GuiValueTypes' is not assignable to type '"text"'.
 *              Type '"dropdown"' is not assignable to type '"text"'.ts(2322)```
 * which is not very helpful.
 */
export type AnyGuiSpec = {
    [T in GuiSpecType]: GuiSpec<T>;
}[GuiSpecType];

export type GuiSpecs = Array<AnyGuiSpec>;
// type TestFunc<T extends GuiSpecType> = GuiSpec<T>['func'];
// const f: TestFunc<'multidropdown'> = (v) => { };
// todo common interface for AddChartDialog & SettingsDialog - from an end-user perspective, not just types
/**
 * This represents widgets used by Add Chart dialog of charts with `extra_controls` specified.
 */
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
// ^^ well-typed as long as we declare the generic type...
// doesn't manage to infer it.
//... we can use g({...}) from `@lib/utils` to get type inference to work. sometimes.

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

export type ChartState<T> = {
    chart: BaseChart<T>;
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
    _popOutChart: (chart: BaseChart) => void;
};
