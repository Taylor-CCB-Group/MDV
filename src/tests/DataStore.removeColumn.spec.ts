import DataStore from "@/datastore/DataStore";
import { describe, expect, test, vi } from "vitest";

function createStore() {
    const tempColumn = {
        field: "temp",
        name: "Temp",
        datatype: "integer",
        data: new Uint8Array([1, 2, 3]),
        buffer: new SharedArrayBuffer(3),
    };

    return Object.assign(Object.create(DataStore.prototype), {
        config: {
            columns: [
                {
                    field: "temp",
                    name: "Temp",
                    datatype: "integer",
                },
                {
                    field: "already_deleted",
                    name: "Already Deleted",
                    datatype: "integer",
                    deleted: true,
                },
            ],
        },
        columnIndex: { temp: tempColumn },
        columns: [tempColumn],
        columnsWithData: ["temp"],
        dirtyColumns: {
            added: {},
            removed: {},
            data_changed: {},
            colors_changed: {},
        },
        dirtyMetadata: new Set(),
        indexes: { temp: { old: 1 } },
        _callListeners: vi.fn(),
    });
}

describe("DataStore soft delete", () => {
    test("excludes soft-deleted columns from the visible runtime store", () => {
        const store = createStore();

        expect(store.columnIndex.temp?.field).toBe("temp");
        expect(store.columnIndex.already_deleted).toBeUndefined();
        expect(store.getAllColumnsMetadata()).toEqual([
            {
                field: "temp",
                name: "Temp",
                datatype: "integer",
            },
            {
                field: "already_deleted",
                name: "Already Deleted",
                datatype: "integer",
                deleted: true,
            },
        ]);
    });

    test("soft-deletes a column through datasource metadata without marking it removed", () => {
        const store = createStore();
        store.dirtyColumns.data_changed.temp = true;
        store.dirtyColumns.colors_changed.temp = true;

        store.softDeleteColumn("temp", true, true);

        expect(store.columnIndex.temp).toBeUndefined();
        expect(store.dirtyColumns.data_changed).toEqual({});
        expect(store.dirtyColumns.colors_changed).toEqual({});
        expect(store.dirtyColumns.removed).toEqual({});
        expect(store.dirtyMetadata.has("columns")).toBe(true);
        expect(store.indexes.temp).toBeUndefined();
        expect(store.getAllColumnsMetadata()).toEqual([
            {
                field: "temp",
                name: "Temp",
                datatype: "integer",
                deleted: true,
            },
            {
                field: "already_deleted",
                name: "Already Deleted",
                datatype: "integer",
                deleted: true,
            },
        ]);
        expect(store._callListeners).toHaveBeenCalledWith("column_removed", "temp");
    });

    test("persists soft-delete metadata even when the column metadata entry is missing", () => {
        const store = createStore();
        store.config.columns = store.config.columns.filter(
            (column: { field: string }) => column.field !== "temp",
        );

        store.softDeleteColumn("temp", true, true);

        expect(store.columnIndex.temp).toBeUndefined();
        expect(store.dirtyMetadata.has("columns")).toBe(true);
        expect(store.getAllColumnsMetadata()).toEqual([
            {
                field: "already_deleted",
                name: "Already Deleted",
                datatype: "integer",
                deleted: true,
            },
            {
                field: "temp",
                name: "Temp",
                datatype: "integer",
                deleted: true,
            },
        ]);
    });

    test("does not add invalid fallback metadata when both runtime and metadata are missing", () => {
        const store = createStore();
        delete store.columnIndex.temp;
        store.columns = [];
        store.columnsWithData = [];
        store.config.columns = store.config.columns.filter(
            (column: { field: string }) => column.field !== "temp",
        );

        store.softDeleteColumn("temp", true, true);

        expect(store.getAllColumnsMetadata()).toEqual([
            {
                field: "already_deleted",
                name: "Already Deleted",
                datatype: "integer",
                deleted: true,
            },
        ]);
        expect(store.dirtyMetadata.has("columns")).toBe(false);
    });
});
