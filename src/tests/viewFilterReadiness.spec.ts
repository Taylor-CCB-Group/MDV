import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import { allSettledWithin, linkCycleMembers, StoredFilters, ViewFilterReadiness } from "@/charts/viewFilterReadiness";

/** Whether `promise` has resolved, without waiting for it beyond the promise reactions already due. */
async function isResolved(promise: Promise<unknown>) {
    let resolved = false;
    promise.then(() => {
        resolved = true;
    });
    for (let i = 0; i < 5; i++) await Promise.resolve();
    return resolved;
}

beforeEach(() => {
    vi.useFakeTimers();
    vi.spyOn(console, "warn").mockImplementation(() => {});
});
afterEach(() => {
    vi.useRealTimers();
    vi.restoreAllMocks();
});

describe("linkCycleMembers", () => {
    test("no dataSource is in a cycle when links only point one way", () => {
        expect(linkCycleMembers({ cells: ["genes"], grid: ["genes"], genes: [] })).toEqual(new Set());
    });

    test("finds the dataSources on a cycle, but not those that only lead into it", () => {
        const members = linkCycleMembers({ a: ["b"], b: ["c"], c: ["b"], d: ["a"] });
        expect(members).toEqual(new Set(["b", "c"]));
    });

    test("a dataSource linked to itself is a cycle", () => {
        expect(linkCycleMembers({ genes: ["genes"], cells: ["genes"] })).toEqual(new Set(["genes"]));
    });
});

describe("allSettledWithin", () => {
    test("resolves to true once everything has settled, whether fulfilled or rejected", async () => {
        const landed = Promise.withResolvers<boolean>();
        const settled = allSettledWithin([landed.promise, Promise.reject(new Error("worker died"))], 5000);
        expect(await isResolved(settled)).toBe(false);
        landed.resolve(true);
        await expect(settled).resolves.toBe(true);
    });

    test("resolves to false at the time limit if something never settles", async () => {
        // e.g. a filter being evaluated on a worker that neither answers nor raises an error
        const settled = allSettledWithin([new Promise(() => {})], 5000);
        await vi.advanceTimersByTimeAsync(4999);
        expect(await isResolved(settled)).toBe(false);
        await vi.advanceTimersByTimeAsync(1);
        await expect(settled).resolves.toBe(false);
    });

    test("leaves no timer running once everything has settled", async () => {
        await allSettledWithin([Promise.resolve()], 5000);
        expect(vi.getTimerCount()).toBe(0);
    });
});

describe("StoredFilters", () => {
    test("with nothing expected, they count as applied", async () => {
        expect(await isResolved(new StoredFilters().whenApplied())).toBe(true);
    });

    test("waits for every expected filter to be reported", async () => {
        const filters = new StoredFilters();
        filters.expect("name");
        filters.expect("score");
        const applied = filters.whenApplied();
        filters.applied("name");
        expect(await isResolved(applied)).toBe(false);
        filters.applied("score");
        expect(await isResolved(applied)).toBe(true);
    });

    test("waits for an asynchronous filter to land", async () => {
        const filters = new StoredFilters();
        filters.expect("score");
        const landed = Promise.withResolvers<boolean>();
        filters.applied("score", landed.promise);
        const applied = filters.whenApplied();
        expect(await isResolved(applied)).toBe(false);
        landed.resolve(true);
        expect(await isResolved(applied)).toBe(true);
    });

    test("a filter that fails to land does not hold things up", async () => {
        const filters = new StoredFilters();
        filters.expect("score");
        filters.applied("score", Promise.reject(new Error("worker died")));
        expect(await isResolved(filters.whenApplied())).toBe(true);
    });

    test("reports for filters that were not expected, or repeated reports, are ignored", async () => {
        const filters = new StoredFilters();
        filters.expect("name");
        filters.applied("other");
        filters.applied("other");
        const applied = filters.whenApplied();
        expect(await isResolved(applied)).toBe(false);
        filters.applied("name");
        filters.applied("name");
        expect(await isResolved(applied)).toBe(true);
    });

    test("gives up waiting after a time limit", async () => {
        const filters = new StoredFilters();
        filters.expect("name");
        const applied = filters.whenApplied(5000);
        await vi.advanceTimersByTimeAsync(4999);
        expect(await isResolved(applied)).toBe(false);
        await vi.advanceTimersByTimeAsync(1);
        expect(await isResolved(applied)).toBe(true);
        expect(console.warn).toHaveBeenCalled();
    });
});

describe("ViewFilterReadiness", () => {
    const dotPlot = { type: "dot_plot" };
    const genesDialog = { type: "selection_dialog" };
    const geneNetwork = { type: "geneNetwork" };

    test("a dataSource with no charts in the view is ready", async () => {
        const readiness = new ViewFilterReadiness({ cells: [dotPlot] });
        expect(await isResolved(readiness.whenReady("genes"))).toBe(true);
    });

    test("a dataSource is ready once all of its charts in the view have settled", async () => {
        const readiness = new ViewFilterReadiness({ cells: [dotPlot], genes: [genesDialog, geneNetwork] });
        const ready = readiness.whenReady("genes");
        readiness.settle("genes", geneNetwork);
        expect(await isResolved(ready)).toBe(false);
        expect(readiness.isReady("genes")).toBe(false);
        readiness.settle("genes", genesDialog);
        expect(readiness.isReady("genes")).toBe(true);
        expect(await isResolved(ready)).toBe(true);
        // ... and stays ready
        expect(await isResolved(readiness.whenReady("genes"))).toBe(true);
    });

    test("charts on other dataSources, and charts that are not part of the view, do not count", async () => {
        const readiness = new ViewFilterReadiness({ cells: [dotPlot], genes: [genesDialog] });
        const ready = readiness.whenReady("genes");
        readiness.settle("cells", dotPlot);
        readiness.settle("genes", { type: "added later by the user" });
        expect(await isResolved(ready)).toBe(false);
    });

    test("a dataSource on a link cycle is not waited for", async () => {
        const readiness = new ViewFilterReadiness({ a: [dotPlot], b: [genesDialog] }, new Set(["a", "b"]));
        expect(await isResolved(readiness.whenReady("b"))).toBe(true);
    });

    test("stops waiting after a time limit, so one stuck chart cannot hold up the rest", async () => {
        const readiness = new ViewFilterReadiness({ genes: [genesDialog] }, new Set(), 10000);
        const ready = readiness.whenReady("genes");
        await vi.advanceTimersByTimeAsync(9999);
        expect(await isResolved(ready)).toBe(false);
        await vi.advanceTimersByTimeAsync(1);
        expect(await isResolved(ready)).toBe(true);
        expect(console.warn).toHaveBeenCalled();
    });

    test("a load that has been abandoned releases whoever was waiting, without a warning", async () => {
        const readiness = new ViewFilterReadiness({ genes: [genesDialog] }, new Set(), 10000);
        const ready = readiness.whenReady("genes");
        readiness.abandon();
        expect(await isResolved(ready)).toBe(true);
        await vi.advanceTimersByTimeAsync(20000);
        expect(console.warn).not.toHaveBeenCalled();
    });

    test("the time limit is not reported when everything settled in time", async () => {
        const readiness = new ViewFilterReadiness({ genes: [genesDialog] }, new Set(), 10000);
        readiness.settle("genes", genesDialog);
        await vi.advanceTimersByTimeAsync(20000);
        expect(console.warn).not.toHaveBeenCalled();
    });
});
