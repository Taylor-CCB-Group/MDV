import { DataModel } from "@/table/DataModel";
import { describe, expect, test, vi } from "vitest";

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

        expect(Array.from(dataStore.columnIndex.label.data)).toEqual([3, 1, 3]);
        expect(dataStore.columnIndex.label.values).toEqual(["", "A", "B", "Filled"]);
        expect(dataChanged).toHaveBeenCalledWith(["label"]);
    });

    test("removes a column through the datastore", () => {
        const removeColumn = vi.fn();
        const dataStore = {
            size: 1,
            columnIndex: {},
            removeColumn,
        } as any;

        const dataModel = new DataModel(dataStore, { autoupdate: false });
        dataModel.removeColumn("label");

        expect(removeColumn).toHaveBeenCalledWith("label", true, true);
    });
});
