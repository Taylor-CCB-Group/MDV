import DataStore from "@/datastore/DataStore";
import { describe, expect, test, vi } from "vitest";
import { type ColumnSyncTarget, type ServerColumn, type ServerDataSource, syncJobColumns } from "./jobColumns";

// DataStore starts a filtering worker when it is built; jsdom has no Worker and these tests never filter
vi.stubGlobal(
    "Worker",
    class {
        postMessage() {}
        terminate() {}
    },
);

const NAME = { field: "name", name: "Name", datatype: "text", values: ["a", "b", "c"] };

function umapColumn(jobId: string, minMax: [number, number]) {
    return {
        field: "UMAP_1",
        name: "UMAP_1",
        datatype: "double",
        minMax,
        provenance: { kind: "job", job_id: jobId, tool_id: "umap", content_hash: "abc" },
    };
}

/** A ChartManager stand-in whose loadColumnSet writes `serverData` into the store, as the real loader does. */
function makeTarget(store: DataStore, serverData: Record<string, number[]> = {}) {
    const loadColumnSet = vi.fn((columns: string[], _dataSource: string, callback: () => void) => {
        for (const c of columns) store.setColumnData(c, serverData[c]);
        callback();
    });
    const target: ColumnSyncTarget = { dataSources: [{ name: "cells", dataStore: store }], loadColumnSet };
    return { target, loadColumnSet };
}

function makeStore(...columns: object[]) {
    const store = new DataStore(3);
    for (const c of columns) store.addColumn(c);
    return store;
}

function column(store: DataStore, field: string) {
    const c = store.columnIndex[field];
    if (!c) throw new Error(`no column ${field}`);
    return c;
}

function server(...columns: ServerColumn[]): ServerDataSource[] {
    return [{ name: "cells", columns: [NAME, ...columns] }];
}

describe("syncJobColumns", () => {
    test("adds a job's new column to the store without loading its data", () => {
        const store = makeStore(NAME);
        const { target, loadColumnSet } = makeTarget(store);

        syncJobColumns(target, server(umapColumn("job1", [0, 1])));

        expect(store.columnIndex.UMAP_1).toMatchObject({ datatype: "double", minMax: [0, 1] });
        expect(store.getAllColumnsMetadata().map((c) => c.field)).toEqual(["name", "UMAP_1"]);
        expect(loadColumnSet).not.toHaveBeenCalled();
    });

    test("reloads a re-run column whose data is loaded and tells charts it changed", () => {
        const store = makeStore(NAME, umapColumn("job1", [0, 1]));
        store.setColumnData("UMAP_1", [0, 0.5, 1]);
        const umap = column(store, "UMAP_1");
        const listener = vi.fn();
        store.addListener("test", listener);
        const { target, loadColumnSet } = makeTarget(store, { UMAP_1: [5, 7, 9] });

        syncJobColumns(target, server(umapColumn("job2", [5, 9])));

        expect(loadColumnSet).toHaveBeenCalledWith(["UMAP_1"], "cells", expect.any(Function));
        // charts holding the column object see the new values
        expect(column(store, "UMAP_1")).toBe(umap);
        expect(umap.data).toEqual(new Float32Array([5, 7, 9]));
        expect(umap).toMatchObject({ minMax: [5, 9], provenance: { job_id: "job2" } });
        expect(listener).toHaveBeenCalledWith("data_changed", expect.objectContaining({ columns: ["UMAP_1"] }));
        expect(store.dirtyColumns.data_changed).not.toHaveProperty("UMAP_1");
    });

    test("updates a re-run column whose data is not loaded without loading it", () => {
        const store = makeStore(NAME, umapColumn("job1", [0, 1]));
        const { target, loadColumnSet } = makeTarget(store);

        syncJobColumns(target, server(umapColumn("job2", [5, 9])));

        expect(store.columnIndex.UMAP_1).toMatchObject({ minMax: [5, 9], provenance: { job_id: "job2" } });
        expect(store.getAllColumnsMetadata().find((c) => c.field === "UMAP_1")).toMatchObject({ minMax: [5, 9] });
        expect(loadColumnSet).not.toHaveBeenCalled();
    });

    test("leaves a column from the same job alone", () => {
        const store = makeStore(NAME, umapColumn("job1", [0, 1]));
        store.setColumnData("UMAP_1", [0, 0.5, 1]);
        const { target, loadColumnSet } = makeTarget(store);

        syncJobColumns(target, server(umapColumn("job1", [0, 1])));

        expect(loadColumnSet).not.toHaveBeenCalled();
        expect(column(store, "UMAP_1").data).toEqual(new Float32Array([0, 0.5, 1]));
    });
});
