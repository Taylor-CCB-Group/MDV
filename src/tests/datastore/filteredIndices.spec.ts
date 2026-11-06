import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import type DataStore from "@/datastore/DataStore";
import { assertConsistent, createStore, type InProcessFilteredIndexWorker, stubWorker } from "./dimensionHarness";

/**
 * `DataStore.getFilteredIndices()`: the request/reply protocol with its worker (running the real worker handler),
 * and what a request resolves to when the filter changes before the worker has got to it.
 */

const ALL = [0, 1, 2, 3, 4, 5, 6, 7];
let ds: DataStore;
let worker: InProcessFilteredIndexWorker;
let dim: ReturnType<DataStore["getDimension"]>;

/** Changes the store's one filter to the rows whose `x` (their own index) is in [min, max] - synchronous at this size. */
function filterTo(min: number, max: number) {
    dim.filter("filterRange", ["x"], { min, max });
}

const indices = async (request: Promise<Uint32Array>) => Array.from(await request);

beforeEach(() => {
    const { filteredIndexWorkers } = stubWorker({ filteredIndex: true });
    vi.spyOn(console, "log").mockImplementation(() => {});
    ds = createStore(8, [{ field: "x", datatype: "double", data: ALL }]);
    worker = filteredIndexWorkers[0];
    dim = ds.getDimension("range_dimension");
});
afterEach(() => {
    assertConsistent(ds);
    vi.restoreAllMocks();
    vi.unstubAllGlobals();
});

describe("DataStore.getFilteredIndices", () => {
    test("resolves with the rows that pass the filter", async () => {
        filterTo(2, 4);
        const request = ds.getFilteredIndices();
        worker.runNext();
        expect(await indices(request)).toEqual([2, 3, 4]);
        expect(ds._filteredIndexPending.size).toBe(0);
    });

    test("calls made while the filter is unchanged share one request", async () => {
        const first = ds.getFilteredIndices();
        const second = ds.getFilteredIndices();
        expect(worker.posted).toBe(1);
        worker.runNext();
        expect(await indices(first)).toEqual(ALL);
        expect(await indices(second)).toEqual(ALL);
        // ... and so does a call made after it has resolved
        expect(await indices(ds.getFilteredIndices())).toEqual(ALL);
        expect(worker.posted).toBe(1);
    });

    test("overlapping requests each settle, with the indices for the current filter", async () => {
        filterTo(2, 4);
        const first = ds.getFilteredIndices();
        filterTo(5, 5); // clears the cache, so the next call is a second request while the first is in flight
        const second = ds.getFilteredIndices();
        expect(worker.posted).toBe(2);

        worker.runNext();
        worker.runNext();

        expect(await indices(second)).toEqual([5]);
        // the first was sized for three rows, and scanned when only row 5 passed: not [5, 0, 0]
        expect(await indices(first)).toEqual([5]);
        expect(ds._filteredIndexPending.size).toBe(0);
    });

    test("a request overtaken by a filter change that nobody has asked about yet is asked again", async () => {
        filterTo(2, 4);
        const request = ds.getFilteredIndices();
        filterTo(5, 6);

        worker.runNext(); // the reply to the overtaken request...
        await Promise.resolve();
        expect(worker.posted).toBe(2); // ... leads to a request for the filter as it is now
        worker.runNext();

        expect(await indices(request)).toEqual([5, 6]);
        expect(ds._filteredIndexPending.size).toBe(0);
    });

    test("if the worker fails, requests in flight and later ones are answered on the main thread", async () => {
        vi.spyOn(console, "error").mockImplementation(() => {});
        filterTo(2, 4);
        const inFlight = ds.getFilteredIndices();

        worker.fail();

        expect(await indices(inFlight)).toEqual([2, 3, 4]);
        expect(ds._filteredIndexPending.size).toBe(0);
        expect(worker.terminated).toBe(true);

        filterTo(6, 7);
        expect(await indices(ds.getFilteredIndices())).toEqual([6, 7]);
        expect(worker.posted).toBe(1);
    });

    test("if posting to the worker throws, it is dropped: that request, ones in flight and later ones are answered on the main thread", async () => {
        vi.spyOn(console, "error").mockImplementation(() => {});
        filterTo(2, 4);
        const inFlight = ds.getFilteredIndices(); // posted, but the worker has not got to it
        filterTo(5, 6);
        const post = vi.spyOn(worker, "postMessage").mockImplementation(() => {
            throw new Error("could not be cloned");
        });

        const failedPost = ds.getFilteredIndices();

        expect(await indices(failedPost)).toEqual([5, 6]);
        expect(await indices(inFlight)).toEqual([5, 6]);
        expect(ds._filteredIndexPending.size).toBe(0);
        expect(worker.terminated).toBe(true);

        filterTo(7, 7);
        expect(await indices(ds.getFilteredIndices())).toEqual([7]);
        expect(post).toHaveBeenCalledOnce();
    });
});
