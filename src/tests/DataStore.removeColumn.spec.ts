import DataStore from "@/datastore/DataStore";
import { describe, expect, test, vi } from "vitest";

describe("DataStore.removeColumn", () => {
    test("clears stale dirty markers for a removed column", () => {
        const column = {
            field: "temp",
            data: new Uint8Array([1, 2, 3]),
            buffer: new SharedArrayBuffer(3),
        };
        const store = {
            columnIndex: { temp: column },
            columns: [column],
            columnsWithData: ["temp"],
            dirtyColumns: {
                added: {},
                removed: {},
                data_changed: { temp: true },
                colors_changed: { temp: true },
            },
            _callListeners: vi.fn(),
        };

        DataStore.prototype.removeColumn.call(store, "temp", true, true);

        expect(store.columnIndex.temp).toBeUndefined();
        expect(store.dirtyColumns.data_changed).toEqual({});
        expect(store.dirtyColumns.colors_changed).toEqual({});
        expect(store.dirtyColumns.removed).toEqual({ temp: true });
        expect(store._callListeners).toHaveBeenCalledWith("column_removed", "temp");
    });

    test("discardPendingAddedColumns removes unsaved added columns without marking them removed", () => {
        const fresh = {
            field: "fresh",
            data: new Uint8Array([1, 2, 3]),
            buffer: new SharedArrayBuffer(3),
        };
        const store = {
            columnIndex: { fresh },
            columns: [fresh],
            columnsWithData: ["fresh"],
            dirtyColumns: {
                added: { fresh: true },
                removed: {},
                data_changed: { fresh: true },
                colors_changed: { fresh: true },
            },
            _callListeners: vi.fn(),
            removeColumn: DataStore.prototype.removeColumn,
        };

        DataStore.prototype.discardPendingAddedColumns.call(store);

        expect(store.columnIndex.fresh).toBeUndefined();
        expect(store.columns).toEqual([]);
        expect(store.columnsWithData).toEqual([]);
        expect(store.dirtyColumns.added).toEqual({});
        expect(store.dirtyColumns.removed).toEqual({});
        expect(store.dirtyColumns.data_changed).toEqual({});
        expect(store.dirtyColumns.colors_changed).toEqual({});
        expect(store._callListeners).not.toHaveBeenCalled();
    });
});
