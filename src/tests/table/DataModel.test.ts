import { DataModel } from "@/table/DataModel";
import { describe, expect, test, vi } from "vitest";

describe("DataModel autoupdate filter listener", () => {
    test("two models on the same store each receive filtered updates", () => {
        const filterArray = new Uint8Array(3);
        filterArray.fill(0);
        const listeners: Record<string, (type: string) => void> = {};
        const dataStore = {
            size: 3,
            filterSize: 3,
            filterArray,
            addListener(id: string, fn: (type: string) => void) {
                listeners[id] = fn;
            },
            removeListener(id: string) {
                delete listeners[id];
            },
        } as any;

        const first = new DataModel(dataStore);
        const second = new DataModel(dataStore);
        first.updateModel();
        second.updateModel();
        expect(Array.from(first.data)).toEqual([0, 1, 2]);
        expect(Array.from(second.data)).toEqual([0, 1, 2]);

        filterArray.set([0, 0, 1]);
        dataStore.filterSize = 2;
        for (const fn of Object.values(listeners)) {
            fn("filtered");
        }
        expect(Array.from(first.data)).toEqual([0, 1]);
        expect(Array.from(second.data)).toEqual([0, 1]);
    });

    test("dispose removes the DataStore listener", () => {
        const filterArray = new Uint8Array(1);
        const listeners: Record<string, (type: string) => void> = {};
        const dataStore = {
            size: 1,
            filterSize: 1,
            filterArray,
            addListener(id: string, fn: (type: string) => void) {
                listeners[id] = fn;
            },
            removeListener(id: string) {
                delete listeners[id];
            },
        } as any;

        const model = new DataModel(dataStore);
        expect(Object.keys(listeners).length).toBe(1);
        model.dispose();
        expect(Object.keys(listeners).length).toBe(0);
        model.dispose();
        expect(Object.keys(listeners).length).toBe(0);
    });
});

describe("DataModel.createColumn", () => {
    test("creates an empty double column using NaN values", () => {
        const addColumn = vi.fn();
        const dataStore = {
            size: 3,
            addColumn,
            columnIndex: {},
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });
        dataModel.createColumn({
            name: "score",
            datatype: "double",
        });

        expect(addColumn).toHaveBeenCalledTimes(1);
        const [column, data, dirty] = addColumn.mock.calls[0];
        expect(column).toMatchObject({
            name: "score",
            field: "score",
            datatype: "double",
            editable: true,
        });
        expect(data).toHaveLength(3);
        expect(data.every((value: number) => Number.isNaN(value))).toBe(true);
        expect(dirty).toBe(true);
    });

    test("creates an empty multitext column with default capacity and delimiter", () => {
        const addColumn = vi.fn();
        const dataStore = {
            size: 2,
            addColumn,
            columnIndex: {},
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });
        dataModel.createColumn({
            name: "tags",
            datatype: "multitext",
        });

        const [column, data] = addColumn.mock.calls[0];
        expect(column).toMatchObject({
            name: "tags",
            field: "tags",
            datatype: "multitext",
            editable: true,
            values: [],
            stringLength: 24,
            delimiter: ",",
        });
        expect(data).toBeInstanceOf(SharedArrayBuffer);
        expect(new Uint16Array(data).every((value) => value === 65535)).toBe(true);
    });

    test("creates an empty unique column with the default string length", () => {
        const addColumn = vi.fn();
        const dataStore = {
            size: 2,
            addColumn,
            columnIndex: {},
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });
        dataModel.createColumn({
            name: "identifier",
            datatype: "unique",
        });

        const [column, data] = addColumn.mock.calls[0];
        expect(column).toMatchObject({
            name: "identifier",
            field: "identifier",
            datatype: "unique",
            editable: true,
            stringLength: 64,
        });
        expect(data).toBeInstanceOf(SharedArrayBuffer);
        expect(data.byteLength).toBe(2 * 64);
    });

    test("clones a unique column preserving datatype metadata and raw data", () => {
        const addColumn = vi.fn();
        const sourceBuffer = new SharedArrayBuffer(6);
        const sourceData = new Uint8Array(sourceBuffer);
        sourceData.set([65, 66, 0, 67, 68, 0]);
        const dataStore = {
            size: 2,
            addColumn,
            columnIndex: {
                source_id: {
                    field: "source_id",
                    name: "Source ID",
                    datatype: "unique",
                    editable: false,
                    stringLength: 3,
                    data: sourceData,
                },
            },
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });
        dataModel.createColumn({
            name: "copied_id",
            cloneColumn: "source_id",
        });

        const [column, data] = addColumn.mock.calls[0];
        expect(column).toMatchObject({
            name: "copied_id",
            field: "copied_id",
            datatype: "unique",
            editable: true,
            stringLength: 3,
        });
        expect(data).toBeInstanceOf(SharedArrayBuffer);
        expect(data).not.toBe(sourceBuffer);
        expect(Array.from(new Uint8Array(data))).toEqual([65, 66, 0, 67, 68, 0]);
    });

    test("creates a compound text column by concatenating source columns", () => {
        const addColumn = vi.fn();
        const dataStore = {
            size: 3,
            addColumn,
            columnIndex: {
                surname: {
                    field: "surname",
                    datatype: "text",
                    data: new Uint8Array([0, 1, 2]),
                    values: ["Jones", "Smith", ""],
                },
                name: {
                    field: "name",
                    datatype: "text",
                    data: new Uint8Array([0, 1, 2]),
                    values: ["Alice", "Bob", ""],
                },
            },
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });
        dataModel.createColumn({
            name: "compound",
            mode: "compound",
            datatype: "text",
            sourceColumns: ["name", "surname"],
            delimiter: "_",
        });

        const [column, data] = addColumn.mock.calls[0];
        expect(column).toMatchObject({
            name: "compound",
            field: "compound",
            datatype: "text",
            editable: true,
        });
        expect(data).toEqual(["Alice_Jones", "Bob_Smith", "_"]);
    });

    test("creates a compound text16 column", () => {
        const addColumn = vi.fn();
        const dataStore = {
            size: 2,
            addColumn,
            columnIndex: {
                a: {
                    field: "a",
                    datatype: "text",
                    data: new Uint8Array([0, 1]),
                    values: ["x", "y"],
                },
                b: {
                    field: "b",
                    datatype: "text",
                    data: new Uint8Array([0, 1]),
                    values: ["u", "v"],
                },
            },
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });
        dataModel.createColumn({
            name: "compound16",
            mode: "compound",
            datatype: "text16",
            sourceColumns: ["a", "b"],
            delimiter: "-",
        });

        const [column, data] = addColumn.mock.calls[0];
        expect(column).toMatchObject({
            name: "compound16",
            field: "compound16",
            datatype: "text16",
            editable: true,
        });
        expect(data).toEqual(["x-u", "y-v"]);
    });

    test("auto-upgrades compound output to text16 when requested text exceeds 256 distinct values", () => {
        const addColumn = vi.fn();
        const size = 300;
        const sourceValues = Array.from({ length: size }, (_, i) => `value_${i}`);
        const sourceData = new Uint16Array(size);
        for (let i = 0; i < size; i++) {
            sourceData[i] = i;
        }
        const dataStore = {
            size,
            addColumn,
            columnIndex: {
                a: {
                    field: "a",
                    datatype: "text16",
                    data: sourceData,
                    values: sourceValues,
                },
                b: {
                    field: "b",
                    datatype: "text",
                    data: new Uint8Array(size),
                    values: [""],
                },
            },
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });
        dataModel.createColumn({
            name: "compound_upgrade",
            mode: "compound",
            datatype: "text",
            sourceColumns: ["a", "b"],
            delimiter: "",
        });

        const [column] = addColumn.mock.calls[0];
        expect(column).toMatchObject({
            datatype: "text16",
        });
    });

    test("throws when compound distinct values exceed text16 capacity", () => {
        const addColumn = vi.fn();
        const size = 65537;
        const sourceValues = Array.from({ length: 65536 }, (_, i) => `value_${i}`);
        const sourceData = new Uint16Array(size);
        for (let i = 0; i < 65536; i++) {
            sourceData[i] = i;
        }
        sourceData[65536] = 0;
        const suffixData = new Uint8Array(size);
        suffixData.fill(0);
        suffixData[65536] = 1;
        const dataStore = {
            size,
            addColumn,
            columnIndex: {
                a: {
                    field: "a",
                    datatype: "text16",
                    data: sourceData,
                    values: sourceValues,
                    getValue: (i: number) => sourceValues[i],
                },
                b: {
                    field: "b",
                    datatype: "text",
                    data: suffixData,
                    values: ["x", "y"],
                    getValue: (i: number) => (i === 65536 ? "y" : "x"),
                },
            },
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });

        expect(() =>
            dataModel.createColumn({
                name: "too_many_distinct",
                mode: "compound",
                datatype: "text16",
                sourceColumns: ["a", "b"],
                delimiter: "",
            }),
        ).toThrow("exceeding text16 limit (65536)");
        expect(addColumn).not.toHaveBeenCalled();
    });
});

describe("DataModel bulk edit helpers", () => {
    test("fills only empty cells in a text column", () => {
        const dataChanged = vi.fn();
        const dataStore = {
            size: 3,
            columnIndex: {
                label: {
                    field: "label",
                    name: "Label",
                    datatype: "text",
                    editable: true,
                    values: ["", "A", "B"],
                    data: new Uint8Array([0, 1, 0]),
                },
            },
            dataChanged,
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });
        dataModel.fillColumn("label", "Filled", [0, 1, 2], true);

        expect(Array.from(dataStore.columnIndex.label.data)).toEqual([1, 0, 1]);
        expect(dataStore.columnIndex.label.values).toEqual(["A", "Filled"]);
        expect(dataChanged).toHaveBeenCalledWith(["label"]);
    });

    test("soft-deletes a column through the datastore", () => {
        const softDeleteColumn = vi.fn();
        const dataStore = {
            size: 1,
            columnIndex: {},
            softDeleteColumn,
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });
        dataModel.removeColumn("label");

        expect(softDeleteColumn).toHaveBeenCalledWith("label", true, true);
    });

    test("throws when softDeleteColumn is unavailable", () => {
        const dataStore = {
            size: 1,
            columnIndex: {},
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });

        expect(() => dataModel.removeColumn("label")).toThrow(
            "DataStore.softDeleteColumn is required for column deletion",
        );
    });

    test("throws for invalid row indices during bulk fill", () => {
        const dataChanged = vi.fn();
        const dataStore = {
            size: 3,
            columnIndex: {
                label: {
                    field: "label",
                    name: "Label",
                    datatype: "text",
                    editable: true,
                    values: [""],
                    data: new Uint8Array([0, 0, 0]),
                },
            },
            dataChanged,
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });

        expect(() => dataModel.fillColumn("label", "Filled", [0, 3], false)).toThrow(
            "Invalid row index 3 for column label",
        );
        expect(dataChanged).not.toHaveBeenCalled();
    });
});
