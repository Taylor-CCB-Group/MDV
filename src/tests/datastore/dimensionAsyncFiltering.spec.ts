import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import type DataStore from "@/datastore/DataStore";
import Dimension from "@/datastore/Dimension";
import { assertConsistent, createStore, passingRows, rng, stubWorker } from "./dimensionHarness";

/**
 * Async filter results: a filter method may hand `commitFilter` a promise (e.g. a mask computed
 * by a worker or GPU kernel). These tests drive resolution order by hand.
 */

type Request = {
    keep: number[];
    resolve: (mask: Uint8Array) => void;
    reject: (e: unknown) => void;
    applied: Promise<boolean>;
};

/** Stands in for a worker/GPU-backed dimension: results arrive when the test resolves them. */
class DeferredDimension extends Dimension {
    requests: Request[] = [];

    filterDeferred(args: { keep: number[] }) {
        const { promise, resolve, reject } = Promise.withResolvers<Uint8Array>();
        const applied = this.commitFilter(promise);
        if (typeof applied === "boolean") throw new Error("expected an async commit");
        this.requests.push({ keep: args.keep, resolve, reject, applied });
        return applied;
    }

    /** Resolves a request with the mask for its own arguments. */
    deliver(request = this.requests[this.requests.length - 1]) {
        const mask = new Uint8Array(this.parent.size);
        for (const i of request.keep) mask[i] = 1;
        request.resolve(mask);
        return request.applied;
    }
}

// row:      0    1    2    3    4    5    6    7
const D = ["p", "p", "q", "q", "p", "q", "p", "q"];
const ALL = [0, 1, 2, 3, 4, 5, 6, 7];

let ds: DataStore;
let events: [string, unknown][];

function deferredDimension() {
    const dim = new DeferredDimension(ds);
    ds.dimensions.push(dim);
    return dim;
}

const settle = () => new Promise((resolve) => setTimeout(resolve, 0));

beforeEach(() => {
    stubWorker();
    vi.spyOn(console, "log").mockImplementation(() => {});
    ds = createStore(8, [
        { field: "x", datatype: "double", data: ALL },
        { field: "d", datatype: "text", data: D },
    ]);
    events = [];
    ds.addListener("spy", (type: string, data: unknown) => events.push([type, data]));
});
afterEach(() => {
    assertConsistent(ds);
    vi.restoreAllMocks();
    vi.unstubAllGlobals();
});

describe("async filter results", () => {
    test("are applied, and listeners notified, only once they resolve", async () => {
        const dim = deferredDimension();
        const done = dim.filter("filterDeferred", [], { keep: [1, 2] });

        expect(dim.pendingFilter).not.toBeNull();
        expect(passingRows(ds)).toEqual(ALL);
        expect(events).toEqual([]);

        dim.deliver();
        await expect(done).resolves.toBe(true);
        expect(passingRows(ds)).toEqual([1, 2]);
        expect(events).toEqual([["filtered", dim]]);
        expect(dim.pendingFilter).toBeNull();
    });

    test("notify=false applies without notifying", async () => {
        const dim = deferredDimension();
        const done = dim.filter("filterDeferred", [], { keep: [1, 2] }, false);
        dim.deliver();
        await expect(done).resolves.toBe(true);
        expect(passingRows(ds)).toEqual([1, 2]);
        expect(events).toEqual([]);
    });

    test("sync filter() also returns a promise, already applied", async () => {
        const dim = ds.getDimension("range_dimension");
        const done = dim.filter("filterRange", ["x"], { min: 2, max: 3 });
        expect(passingRows(ds)).toEqual([2, 3]);
        expect(dim.pendingFilter).toBeNull();
        await expect(done).resolves.toBe(true);
    });

    test("an older result that lands first is shown until the latest arrives", async () => {
        const dim = deferredDimension();
        const first = dim.filter("filterDeferred", [], { keep: [1, 2] });
        const second = dim.filter("filterDeferred", [], { keep: [5] });

        dim.deliver(dim.requests[0]);
        await expect(first).resolves.toBe(true);
        expect(passingRows(ds)).toEqual([1, 2]);
        expect(dim.pendingFilter).not.toBeNull();

        dim.deliver(dim.requests[1]);
        await expect(second).resolves.toBe(true);
        expect(passingRows(ds)).toEqual([5]);
        expect(events).toEqual([
            ["filtered", dim],
            ["filtered", dim],
        ]);
    });

    test("an older result that lands after a newer one is dropped", async () => {
        const dim = deferredDimension();
        const first = dim.filter("filterDeferred", [], { keep: [1, 2] });
        const second = dim.filter("filterDeferred", [], { keep: [5] });
        dim.deliver(dim.requests[1]);
        dim.deliver(dim.requests[0]);

        await expect(second).resolves.toBe(true);
        await expect(first).resolves.toBe(false);
        expect(passingRows(ds)).toEqual([5]);
        expect(events).toEqual([["filtered", dim]]);
    });

    test("adding rows cancels results computed for the old row count", async () => {
        const dim = deferredDimension();
        const stale = dim.filter("filterDeferred", [], { keep: [1, 2] });
        ds.addDataToStore(
            [
                { column: "x", data: new Float32Array([8, 9]).buffer },
                { column: "d", data: new Uint8Array([0, 0]).buffer },
            ],
            2,
        );
        expect(dim.requests).toHaveLength(2);

        dim.deliver(dim.requests[0]);
        await expect(stale).resolves.toBe(false);
        await dim.deliver(dim.requests[1]);
        expect(passingRows(ds)).toEqual([1, 2]);
    });

    test("a result resolving to null is not applied", async () => {
        const dim = deferredDimension();
        const pending = dim.filter("filterDeferred", [], { keep: [1] });
        dim.requests[0].resolve(null as never);
        await expect(pending).resolves.toBe(false);
        expect(passingRows(ds)).toEqual(ALL);
    });

    test("a sync filter supersedes a pending async one", async () => {
        const dim = deferredDimension();
        const pending = dim.filter("filterDeferred", [], { keep: [1, 2] });
        dim.filter("filterPredicate", [], { predicate: (i: number) => i > 5 });
        dim.deliver(dim.requests[0]);
        await expect(pending).resolves.toBe(false);
        expect(passingRows(ds)).toEqual([6, 7]);
    });

    test("removeFilter drops a pending result", async () => {
        const dim = deferredDimension();
        const pending = dim.filter("filterDeferred", [], { keep: [1, 2] });
        dim.removeFilter();
        expect(dim.pendingFilter).toBeNull();
        dim.deliver();
        await expect(pending).resolves.toBe(false);
        expect(passingRows(ds)).toEqual(ALL);
    });

    test("removeAllFilters drops pending results", async () => {
        const dim = deferredDimension();
        const pending = dim.filter("filterDeferred", [], { keep: [1, 2] });
        ds.removeAllFilters();
        dim.deliver();
        await expect(pending).resolves.toBe(false);
        expect(passingRows(ds)).toEqual(ALL);
    });

    test("destroy drops a pending result without touching the parent", async () => {
        const other = ds.getDimension("range_dimension");
        other.filter("filterRange", ["x"], { min: 0, max: 5 });
        const dim = deferredDimension();
        const pending = dim.filter("filterDeferred", [], { keep: [1, 2] });
        dim.destroy();
        dim.deliver();
        await expect(pending).resolves.toBe(false);
        expect(passingRows(ds)).toEqual([0, 1, 2, 3, 4, 5]);
    });

    test("a short mask is rejected without changing state", async () => {
        const dim = deferredDimension();
        const pending = dim.filter("filterDeferred", [], { keep: [] });
        dim.requests[0].resolve(new Uint8Array(3));
        await expect(pending).rejects.toThrow(/3 rows, expected 8/);
        expect(passingRows(ds)).toEqual(ALL);
        expect(dim.pendingFilter).toBeNull();
    });

    test("a rejected evaluation leaves the previous result in place", async () => {
        const dim = deferredDimension();
        dim.filter("filterDeferred", [], { keep: [1, 2] });
        await dim.deliver();
        const failed = dim.filter("filterDeferred", [], { keep: [3] });
        dim.requests[1].reject(new Error("device lost"));
        await expect(failed).rejects.toThrow("device lost");
        expect(passingRows(ds)).toEqual([1, 2]);
    });

    test("changing the background filter re-requests and notifies once the result lands", async () => {
        const dim = deferredDimension();
        dim.filter("filterDeferred", [], { keep: [0, 1, 2, 3] });
        await dim.deliver();
        events.length = 0;

        const done = dim.setBackgroundFilter("d", "p");
        expect(dim.requests).toHaveLength(2);
        expect(events).toEqual([]);
        assertConsistent(ds);

        dim.deliver();
        await expect(done).resolves.toBe(true);
        // rows outside "p" are background-hidden and, with a local filter active, excluded
        expect(Array.from(dim.filterArray)).toEqual([0, 0, 3, 3, 1, 3, 1, 3]);
        expect(passingRows(ds)).toEqual([0, 1]);
        expect(events).toEqual([["filtered", dim]]);

        const cleared = dim.clearBackgroundFilter();
        dim.deliver();
        await expect(cleared).resolves.toBe(true);
        expect(passingRows(ds)).toEqual([0, 1, 2, 3]);
    });

    test("dataChanged reports async refilters once they land", async () => {
        const dim = deferredDimension();
        dim.filter("filterDeferred", ["x"], { keep: [1] });
        await dim.deliver();
        events.length = 0;

        ds.dataChanged(["x"], false);
        expect(events).toEqual([["data_changed", { columns: ["x"], hasFiltered: false }]]);

        dim.deliver();
        await settle();
        expect(events.map(([type]) => type)).toEqual(["data_changed", "filtered"]);
    });

    test("dataChanged still notifies when one async refilter fails", async () => {
        const error = vi.spyOn(console, "error").mockImplementation(() => {});
        const ok = deferredDimension();
        const failing = deferredDimension();
        ok.filter("filterDeferred", ["x"], { keep: [1] });
        failing.filter("filterDeferred", ["x"], { keep: [2] });
        await Promise.all([ok.deliver(), failing.deliver()]);
        events.length = 0;

        ds.dataChanged(["x"], false);
        failing.requests[1].reject(new Error("device lost"));
        ok.deliver();
        await settle();

        expect(events.map(([type]) => type)).toEqual(["data_changed", "filtered"]);
        expect(error).toHaveBeenCalledWith("async refilter failed", expect.any(Error));
        expect(passingRows(ds)).toEqual([]);
    });
});

/**
 * Random interleavings of async requests, resolutions (in any order), cancellations and
 * background changes. Every resolution must be applied iff it is newer than the applied result
 * and wasn't cancelled; counts must stay consistent at every step; once nothing is pending, each
 * dimension must hold exactly the latest live request's result.
 */
describe("randomized async interleavings", () => {
    type Model = {
        applied: Set<number> | null;
        live: Request | null;
        active: boolean;
        bg: Set<string> | null;
        /** request order, as commitFilter numbers them */
        seq: number;
        appliedSeq: number;
        cancelledThrough: number;
    };

    function expectedState(m: Model, i: number) {
        if (m.bg && !m.bg.has(D[i])) return m.applied ? 3 : 2;
        return m.applied && !m.applied.has(i) ? 1 : 0;
    }

    test.each([1, 2, 3, 4, 5, 6, 7, 8])("seed %i", async (seed) => {
        const rand = rng(seed);
        const dims = [deferredDimension(), deferredDimension(), deferredDimension()];
        const models: Model[] = dims.map(() => ({
            applied: null,
            live: null,
            active: false,
            bg: null,
            seq: 0,
            appliedSeq: 0,
            cancelledThrough: 0,
        }));
        const outstanding = new Set<[DeferredDimension, Model, Request, number]>();

        const track = (i: number) => {
            const dim = dims[i];
            const request = dim.requests[dim.requests.length - 1];
            const m = models[i];
            m.live = request;
            outstanding.add([dim, m, request, ++m.seq]);
        };
        const cancel = (m: Model) => {
            m.live = null;
            m.applied = null;
            m.active = false;
            m.cancelledThrough = m.seq;
        };

        for (let step = 0; step < 80; step++) {
            const roll = rand();
            const i = Math.floor(rand() * dims.length);
            const [dim, m] = [dims[i], models[i]];
            if (roll < 0.3) {
                const keep = ALL.filter(() => rand() < 0.5);
                void dim.filter("filterDeferred", [], { keep });
                m.active = true;
                track(i);
            } else if (roll < 0.6 && outstanding.size) {
                const all = Array.from(outstanding);
                const entry = all[Math.floor(rand() * all.length)];
                outstanding.delete(entry);
                const [d, model, request, seq] = entry;
                const shouldApply = seq > model.appliedSeq && seq > model.cancelledThrough;
                await expect(d.deliver(request), `step ${step} seq ${seq}`).resolves.toBe(shouldApply);
                if (shouldApply) {
                    model.applied = new Set(request.keep);
                    model.appliedSeq = seq;
                }
                if (model.live === request) model.live = null;
            } else if (roll < 0.7) {
                dim.removeFilter();
                cancel(m);
            } else if (roll < 0.74) {
                ds.removeAllFilters();
                models.forEach(cancel);
            } else if (roll < 0.9) {
                const cats = rand() < 0.5 ? ["p"] : ["q"];
                void dim.setBackgroundFilter("d", cats);
                m.bg = new Set(cats);
                if (m.active) track(i);
            } else {
                const hadBg = m.bg !== null;
                void dim.clearBackgroundFilter();
                m.bg = null;
                if (hadBg && m.active) track(i);
            }
            await settle();
            assertConsistent(ds);
            dims.forEach((d, n) => {
                expect(d.pendingFilter === null, `step ${step} dim ${n} pending`).toBe(models[n].live === null);
                if (models[n].live === null) {
                    const expected = ALL.map((row) => expectedState(models[n], row));
                    expect(Array.from(d.filterArray), `step ${step} dim ${n}`).toEqual(expected);
                }
            });
        }
    });
});
