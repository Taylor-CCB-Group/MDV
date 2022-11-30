export type DataSource = { contentDiv: HTMLDivElement };

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
