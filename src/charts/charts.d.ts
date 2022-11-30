// Not necessarily comprehensive set of types, arbitrary coverage of
// - parts of code that intersect with stuff I'm actually using directly in TS
// - general shape of data structures 

export type DataType = 'integer' | 'double' | 'text' | 'unique';

type DataStructureTypes = {
    'integer': Uint32Array;
    'double': Float32Array; //why is it called 'double'???
    'text': Uint8Array;
    'unique': Uint8Array; //not sure about this either.
}
type DataValuesTypes = {
    'integer': undefined;
    'double': undefined;
    'text': string[];
    'unique': string[]; //not totally sure
}

export type DataColumn<T extends DataType> = {
    name: string;
    field: string;
    datatype: T;
    data: DataStructureTypes[T];
    values: DataValuesTypes[T];
}


export type DataStore = {
    size: number;
    filterSize: number;
    columns: DataColumn<any>[];
    //tempted to try to be clever about generics in strings (keyof etc)
    columnGroups: Record<string, {name: string, columns: string[]}>;
    columnIndex: Record<string, DataColumn<any>>;
    columnsWithData: string[];
};
export type DataSource = {
    name: string;
    charts: Chart[];
    dataStore: DataStore;
    contentDiv: HTMLDivElement;
    menuBar: HTMLDivElement;
};

export type GuiSpec<T> = {
    type: string; 
    label: string;
    current_value: T;
    func: (T) => void;
    values?: T[];
}

export type Chart = {
    getDiv: () => HTMLDivElement;
    remove: () => void;
    addMenuIcon: (classes: string, info: string) => HTMLElement;
    setSize: (x?: number, y?: number) => void;
    changeBaseDocument: (doc: Document) => void;
    getSettings: () => GuiSpec<any>[]
};

export type ChartManager = {
    addMenuIcon: (dataSourceName: string, iconClass: string, text: string, func: ()=>void) => void;
};
