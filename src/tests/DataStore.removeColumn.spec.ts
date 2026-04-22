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
});
