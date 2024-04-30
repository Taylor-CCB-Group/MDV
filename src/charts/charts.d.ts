// Not necessarily comprehensive set of types, arbitrary coverage of
// - parts of code that intersect with stuff I'm actually using directly in TS
// - general shape of data structures NB - in many cases this can be inferred from JS, so I may remove these again...
//   partly just using them as a form of documentation / notes-to-self.

import DataStore from "../datastore/DataStore";

export type DataType = 'integer' | 'double' | 'text' | 'unique' | 'multitext' | 'int32';

type DataStructureTypes = {
    'integer': Uint32Array;
    'double': Float32Array; //why is it called 'double'???
    'text': Uint8Array;
    'multitext': Uint16Array;
    'unique': Uint8Array; //not sure about this either.
}
type DataValuesTypes = {
    'integer': undefined;
    'double': undefined;
    'text': string[]; //would be better if this was `Set<string>`? maybe not, want indexOf
    'multitext': string[]; //would be better if this was `Set<string>`? maybe not, want indexOf
    'unique': string[];
}
// even if they're just aliases, these could be useful for documentation / clarity
export type ColumnName = string;
export type DataSourceName = string;
export type FieldName = string;

export type DataColumn<T extends DataType> = {
    name: ColumnName;
    field: FieldName;
    datatype: T;
    data: DataStructureTypes[T];
    values: DataValuesTypes[T];
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
export type DataSource = {
    name: DataSourceName;
    charts: Chart[];
    dataStore: DataStore;
    contentDiv: HTMLDivElement;
    menuBar: HTMLDivElement;
    images?: Record<string, any>;
    regions?: Record<string, any>;
};

export type DropdownMappedValue<T extends string, V extends string> = { [P in T]: string } & { [P in V]: string };
export type DropdownMappedValues<T extends string, V extends string> = [
    Array<DropdownMappedValue<T, V>>, T, V
];
export type DropDownValues = DropdownMappedValues<string, string> | [Array<string>];
export type GuiValueTypes = {
    "dropdown": string;
    "multidropdown": any;
    "check": boolean;
    "text": string;
    "textbox": string;
    "radiobuttons": string;
    "slider": number;
    "spinner": number;
    "button": undefined;
    "doubleslider": [number, number];
    "folder": GuiSpec[];
}
export type GuiSpecType = keyof GuiValueTypes;
export type GuiSpec<T extends GuiSpecType> = {
    type: T; 
    label: string;
    name: string;
    current_value?: GuiValueTypes[T];
    func?: (v: GuiValueTypes[T]) => void;
    values?: T extends ('dropdown' | 'multidropdown') ? DropDownValues : never;
    // choices is only used for radiobuttons, so we should infer if T is radiobuttons, otherwise never
    choices?: T extends 'radiobuttons' ? [string, string][] : never;
    min?: number;
    max?: number;
    step?: number;
    defaultVal?: GuiValueTypes[T];
}
// todo common interface for AddChartDialog & SettingsDialog - from an end-user perspective, not just types
export type ExtraControl<T extends GuiSpecType> = {
    type: T;
    name: string;
    label: string;
    values?: Array<{name: string, value: string}>;
    defaultVal?: GuiValueTypes[T];
}
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
interface DataStore {
    getLoadedColumns: () => FieldName[];
    getColumnName: (col: FieldName) => ColumnName | null;
    getColumnList: () => DataColumn[];
    getFilteredIndices: () => Promise<Uint32Array>;
}


export interface Chart {
    getDiv: () => HTMLElement;
    remove: () => void;
    addMenuIcon: (classes: string, info: string) => HTMLElement;
    setSize: (x?: number, y?: number) => void;
    changeBaseDocument: (doc: Document) => void;
    getSettings: () => GuiSpec[];
    removeLayout?:()=> void;
    config:any;
    dataStore: DataStore;
};    

export type ChartState = {
    chart: Chart;
    win?: Window;
    dataSource: DataSource;
}    

export type DataSourceSpec = {
    name: DataSourceName;
    dataStore: DataStore;
    //links etc TBD
};

export type ChartManager = {
    charts: Record<string, ChartState>;
    dataSources: DataSourceSpec[];
    dsIndex: Record<string, DataSourceSpec>;
    addMenuIcon: (dataSourceName: DataSourceName, iconClass: string, text: string, func: ()=>void) => HTMLElement;
    /** probably not something we really want to use publicly like this */
    _popOutChart: (chart: Chart) => void;
};
