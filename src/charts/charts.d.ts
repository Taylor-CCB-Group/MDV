export type DataStore = {
    size: number, filterSize: number; columns: any[]; //etc
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
