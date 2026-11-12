import { afterEach, describe, expect, test, vi } from "vitest";
import { withLoadedColumns } from "@/charts/chartConfigUtils";
import type DataStore from "@/datastore/DataStore";
import { RowsAsColsQuery } from "@/links/link_utils";

const query = { linkedDsName: "rna", maxItems: 10, type: "RowsAsColsQuery" as const };
const ds = {} as DataStore;

/** Simulates the column data held by the DataStore, keyed by field name. */
function fakeStore(initiallyLoaded: string[]) {
    const loaded = new Set(initiallyLoaded);
    return {
        loaded,
        getMissing: (fields: string[]) => fields.filter((f) => !loaded.has(f)),
        /** `load` as ChartManager would call it: given specs, resolves the fields as they are now */
        loadSpecs: (resolve: (specs: any[]) => string[]) => async (specs: any[]) => {
            for (const f of resolve(specs)) loaded.add(f);
        },
    };
}

/** A query whose fields are whatever `current.value` says - to simulate a filter landing at an awkward moment. */
function liveQuery(current: { value: string }) {
    return {
        get fields() {
            return [current.value];
        },
        async initialize() {},
    };
}

afterEach(() => vi.restoreAllMocks());

describe("withLoadedColumns", () => {
    test("loads the new field when a query's fields change while loading", async () => {
        const gene = { value: "rna_expr|OLD (rna_expr)|1" };
        vi.spyOn(RowsAsColsQuery, "fromSerialized").mockReturnValue(liveQuery(gene) as any);
        const store = fakeStore(["rna:bucket"]);
        const load = vi.fn(async () => {
            store.loaded.add(gene.value);
            // the linked gene filter lands while this load is in flight
            gene.value = "rna_expr|NEW (rna_expr)|96";
        });

        // `use` stands in for the chart constructor, which resolves the query again
        const seen = await withLoadedColumns(
            ds,
            () => ["rna:bucket", query],
            store.getMissing,
            load,
            () => gene.value,
        );

        expect(seen).toBe("rna_expr|NEW (rna_expr)|96");
        expect(store.getMissing([seen])).toEqual([]);
        expect(load).toHaveBeenCalledTimes(2);
    });

    test("covers queries outside `param` (e.g. tooltip.column) that change while loading", async () => {
        const param = { value: "rna_expr|A (rna_expr)|1" };
        const tooltip = { value: "rna_expr|B (rna_expr)|2" };
        // `getSpecs()` returns [param query, tooltip query], deserialised in that order on every attempt
        let calls = 0;
        vi.spyOn(RowsAsColsQuery, "fromSerialized").mockImplementation(
            () => liveQuery(calls++ % 2 === 0 ? param : tooltip) as any,
        );
        const store = fakeStore([]);
        const load = vi.fn(async () => {
            store.loaded.add(param.value);
            store.loaded.add(tooltip.value);
            // only the first load is interrupted by the tooltip's gene changing
            if (tooltip.value === "rna_expr|B (rna_expr)|2") tooltip.value = "rna_expr|C (rna_expr)|3";
        });

        const seen = await withLoadedColumns(
            ds,
            () => [query, query],
            store.getMissing,
            load,
            () => [param.value, tooltip.value],
        );

        expect(seen).toEqual(["rna_expr|A (rna_expr)|1", "rna_expr|C (rna_expr)|3"]);
        expect(store.getMissing(seen)).toEqual([]);
    });

    test("waits for a query to initialise before reading its fields", async () => {
        const gene = { value: "rna_expr|G (rna_expr)|1" };
        let initialised = false;
        const q = {
            get fields() {
                if (!initialised) throw new Error("fields read before initialize()");
                return [gene.value];
            },
            initialize: async () => {
                await Promise.resolve();
                initialised = true;
            },
        };
        vi.spyOn(RowsAsColsQuery, "fromSerialized").mockReturnValue(q as any);
        const store = fakeStore([gene.value]);
        const use = vi.fn(() => "chart");
        await expect(withLoadedColumns(ds, () => [query], store.getMissing, vi.fn(), use)).resolves.toBe("chart");
        expect(use).toHaveBeenCalledOnce();
    });

    test("a change that is already queued when the check passes cannot get in before `use`", async () => {
        const gene = { value: "rna_expr|OLD (rna_expr)|1" };
        vi.spyOn(RowsAsColsQuery, "fromSerialized").mockReturnValue(liveQuery(gene) as any);
        const store = fakeStore([gene.value]);
        const getMissing = (fields: string[]) => {
            // stands in for a promise reaction (e.g. a filter result being applied to the link) that is already
            // queued: it runs before anything that resumes from an `await` on this function's result
            queueMicrotask(() => {
                gene.value = "rna_expr|NEW (rna_expr)|96";
            });
            return store.getMissing(fields);
        };

        const seen = await withLoadedColumns(
            ds,
            () => [query],
            getMissing,
            vi.fn(),
            () => gene.value,
        );

        expect(store.getMissing([seen])).toEqual([]);
        // ... whereas by the time the caller resumes, the query has moved on to a column that isn't loaded
        expect(store.getMissing([gene.value])).toEqual([gene.value]);
    });

    test("does not load anything when everything is already loaded", async () => {
        const store = fakeStore(["a", "b"]);
        const load = vi.fn();
        const use = vi.fn();
        await withLoadedColumns(ds, () => ["a", "b"], store.getMissing, load, use);
        expect(load).not.toHaveBeenCalled();
        expect(use).toHaveBeenCalledOnce();
    });

    test("uses columns that only become available on the last permitted load", async () => {
        const store = fakeStore([]);
        const load = vi.fn(async () => {
            store.loaded.add("a");
        });
        const use = vi.fn(() => "chart");
        await expect(withLoadedColumns(ds, () => ["a"], store.getMissing, load, use, 1)).resolves.toBe("chart");
        expect(load).toHaveBeenCalledOnce();
    });

    test("checks the specs the config has after initialisation, not the ones it had before", async () => {
        // the config is given a different column while we wait for a query to initialise
        let specs: (string | typeof query)[] = [query];
        const gene = { value: "rna_expr|G (rna_expr)|1" };
        const q = {
            get fields() {
                return [gene.value];
            },
            initialize: async () => {
                specs = ["other"];
            },
        };
        vi.spyOn(RowsAsColsQuery, "fromSerialized").mockReturnValue(q as any);
        const store = fakeStore([gene.value]);
        const load = vi.fn(async (toLoad: (string | typeof query)[]) => {
            for (const spec of toLoad) if (typeof spec === "string") store.loaded.add(spec);
        });

        // `use` stands in for the chart constructor, which reads the config as it is by then
        const seen = await withLoadedColumns(
            ds,
            () => specs,
            store.getMissing,
            load,
            () => specs,
        );

        expect(seen).toEqual(["other"]);
        expect(store.getMissing(["other"])).toEqual([]);
    });

    test("gives up if the query never settles", async () => {
        let n = 0;
        const q = {
            get fields() {
                return [`rna_expr|G${n++} (rna_expr)|${n}`];
            },
            async initialize() {},
        };
        vi.spyOn(RowsAsColsQuery, "fromSerialized").mockReturnValue(q as any);
        const store = fakeStore([]);
        const use = vi.fn();
        await expect(withLoadedColumns(ds, () => [query], store.getMissing, vi.fn(), use, 3)).rejects.toThrow(
            /kept changing/,
        );
        expect(use).not.toHaveBeenCalled();
    });
});
