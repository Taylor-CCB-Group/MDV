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

export type DataColumn<T extends DataType> = {
    name: string;
    field: string;
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
    name: string;
    charts: Chart[];
    dataStore: DataStore;
    contentDiv: HTMLDivElement;
    menuBar: HTMLDivElement;
};

export type GuiValueTypes = {
    "dropdown": string;
    "check": boolean;
    "text": string;
    "radiobuttons": any;
    "slider": number;
    "button": undefined; //never?
    "doubleslider": [number, number];
}

export type GuiSpec<T extends keyof GuiValueTypes> = {
    type: T; 
    label: string;
    current_value: GuiValueTypes[T];
    func: (v: GuiValueTypes[T]) => void;
    values?: GuiValueTypes[T][];
}

export type Chart = {
    getDiv: () => HTMLDivElement;
    remove: () => void;
    addMenuIcon: (classes: string, info: string) => HTMLElement;
    setSize: (x?: number, y?: number) => void;
    changeBaseDocument: (doc: Document) => void;
    getSettings: () => GuiSpec<any>[];
    removeLayout:()=> void;
    config:any;
};

export type ChartManager = {
    addMenuIcon: (dataSourceName: string, iconClass: string, text: string, func: ()=>void) => HTMLElement;
};
