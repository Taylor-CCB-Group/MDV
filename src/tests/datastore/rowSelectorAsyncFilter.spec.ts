import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import type DataStore from "@/datastore/DataStore";
import Dimension from "@/datastore/Dimension";
import { createRowSelector } from "@/links/link_utils";
import { assertConsistent, createStore, passingRows, stubWorker } from "./dimensionHarness";

/**
 * `createRowSelector` decides whether a highlighted row is visible by reading `ds.filterArray`.
 * Range filters can now be evaluated on a worker, landing later, so these tests check that reading it is
 * safe against a real DataStore / Dimension whose async results we deliver by hand.
 */

/** Stands in for a worker-backed dimension: a filter's result is applied when the test delivers it. */
class DeferredDimension extends Dimension {
    resolvers: ((mask: Uint8Array) => void)[] = [];
    /** A filter method that keeps only the rows in `keep`, once delivered. */
    filterDeferred(args: { keep: number[] }) {
        const { promise, resolve } = Promise.withResolvers<Uint8Array>();
        this.resolvers.push((mask) => {
            for (const i of args.keep) mask[i] = 1;
            resolve(mask);
        });
        return this.commitFilter(promise);
    }
    /** Resolves a pending filter result (the latest by default), as the worker replying would. */
    deliver(index = this.resolvers.length - 1) {
        this.resolvers[index](new Uint8Array(this.parent.size));
    }
}

const ALL = [0, 1, 2, 3, 4, 5, 6, 7];
let ds: DataStore;

/** the real store's highlight and filter state, with the worker-backed `getFilteredIndices` replaced */
function selectorFor(ds: DataStore, applied: number[][]) {
    const source = {
        getHighlightedData: () => ds.getHighlightedData(),
        get filterArray() {
            return ds.filterArray;
        },
        getFilteredIndices: async () => new Uint32Array(passingRows(ds)),
    };
    return createRowSelector(source, (rows) => applied.push(Array.from(rows)));
}

beforeEach(() => {
    stubWorker();
    vi.spyOn(console, "log").mockImplementation(() => {});
    ds = createStore(8, [{ field: "x", datatype: "double", data: ALL }]);
});
afterEach(() => {
    assertConsistent(ds);
    vi.restoreAllMocks();
    vi.unstubAllGlobals();
});

describe("createRowSelector against async filters", () => {
    test("filterArray already reflects an async filter when 'filtered' is raised", async () => {
        const dim = new DeferredDimension(ds);
        ds.dimensions.push(dim);
        let rowsAtEvent: number[] | undefined;
        ds.addListener("spy", (type: string) => {
            if (type === "filtered") rowsAtEvent = passingRows(ds);
        });
        const done = dim.filter("filterDeferred", [], { keep: [1, 2] });
        dim.deliver();
        await done;
        expect(rowsAtEvent).toEqual([1, 2]);
    });

    test("filterArray does not change before an async filter lands", async () => {
        const dim = new DeferredDimension(ds);
        ds.dimensions.push(dim);
        const done = dim.filter("filterDeferred", [], { keep: [1, 2] });
        await new Promise((r) => setTimeout(r, 0));
        expect(passingRows(ds)).toEqual(ALL);
        dim.deliver();
        await done;
        expect(passingRows(ds)).toEqual([1, 2]);
    });

    test("a highlight is judged against the filter as committed, then re-judged when a filter lands", async () => {
        const dim = new DeferredDimension(ds);
        ds.dimensions.push(dim);
        const applied: number[][] = [];
        const refresh = selectorFor(ds, applied);
        ds.addListener("link", (type: string) => {
            if (type === "filtered") refresh();
        });
        ds.highightedData = [3];

        await refresh();
        expect(applied.at(-1)).toEqual([3]);

        // filter in flight: row 3 still passes the committed filter, so the highlight still stands
        const done = dim.filter("filterDeferred", [], { keep: [1, 2] });
        await refresh();
        expect(applied.at(-1)).toEqual([3]);

        // once it lands, 'filtered' is raised, and row 3 is no longer visible: fall back to the filtered rows
        dim.deliver();
        await done;
        await new Promise((r) => setTimeout(r, 0));
        expect(applied.at(-1)).toEqual([1, 2]);
    });

    // `filterArray` changes when an async result lands, but 'filtered' is raised a few microtask ticks later. A
    // highlight event can arrive inside that window (it is raised from a promise), and if the highlight is filtered
    // out the selector asks for the filtered rows: those must be the rows for the committed filter, not the ones
    // cached for the filter before it.
    test.each([0, 1, 2, 3, 4, 5, 6, 7])(
        "a highlight event %i ticks after an async filter commits still ends consistent",
        async (ticks) => {
            const dim = new DeferredDimension(ds);
            ds.dimensions.push(dim);
            const applied: number[][] = [];
            const source = {
                getHighlightedData: () => ds.getHighlightedData(),
                get filterArray() {
                    return ds.filterArray;
                },
                // same caching behaviour as DataStore.getFilteredIndices, minus the worker
                getFilteredIndices: async () => {
                    if (!ds._filteredIndicesPromise) {
                        const indices = Promise.resolve(new Uint32Array(passingRows(ds)));
                        Object.assign(ds, { _filteredIndicesPromise: indices });
                        return indices;
                    }
                    return ds._filteredIndicesPromise;
                },
            };
            const refresh = createRowSelector(source, (rows) => applied.push(Array.from(rows)));
            ds.addListener("link", (type: string) => {
                if (type === "filtered" || type === "data_highlighted") refresh();
            });
            ds.highightedData = [3];
            await source.getFilteredIndices(); // cache the unfiltered rows, as an earlier consumer would have
            await refresh();
            const done = dim.filter("filterDeferred", [], { keep: [1, 2] });
            dim.deliver();
            for (let i = 0; i < ticks; i++) await Promise.resolve();
            ds._callListeners("data_highlighted", { indexes: [3] });
            await done;
            await new Promise((r) => setTimeout(r, 0));
            expect(applied.at(-1)).toEqual([1, 2]);
            // and the rows cached for the filter before the commit are never exposed after it
            expect(applied.slice(1)).not.toContainEqual(ALL);
        },
    );
});
