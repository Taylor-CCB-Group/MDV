import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import type DataStore from "@/datastore/DataStore";
import RangeDimension from "@/datastore/RangeDimension";
import { evaluateFilterSpec, type FilterSpec, type NumericArray } from "@/datastore/rangeFilterEval";
import { RangeFilterClient } from "@/datastore/rangeFilterClient";
import {
    assertConsistent,
    createStore,
    InProcessFilterWorker,
    passingRows,
    rng,
    settleFilters,
    stubWorker,
} from "./dimensionHarness";

/**
 * RangeDimension's worker path: the pure evaluator must match the predicates RangeDimension used
 * before it existed, and the async path must commit the same states as the sync one.
 */

/** The pre-worker RangeDimension predicates, verbatim apart from resolving columns up front. */
const legacy = {
    range(arr: NumericArray, min: number, max: number) {
        return (i: number) => {
            const v = arr[i];
            return v >= min && v <= max && !Number.isNaN(v);
        };
    },
    square(data1: NumericArray, data2: NumericArray, range1: number[], range2: number[]) {
        return (i: number) => {
            const v1 = data1[i];
            const v2 = data2[i];
            return (
                v1 >= range1[0] &&
                v1 <= range1[1] &&
                v2 >= range2[0] &&
                v2 <= range2[1] &&
                Number.isFinite(v1) &&
                Number.isFinite(v2)
            );
        };
    },
    poly(
        data1: NumericArray,
        data2: NumericArray,
        points: [number, number][],
        isRegionScoped: boolean,
        regionData: NumericArray | null,
        wantedRegionIndex: number,
    ) {
        let minX = Number.MAX_VALUE;
        let minY = Number.MAX_VALUE;
        let maxX = Number.MIN_VALUE;
        let maxY = Number.MIN_VALUE;
        for (const pt of points) {
            minX = Math.min(minX, pt[0]);
            maxX = Math.max(maxX, pt[0]);
            minY = Math.min(minY, pt[1]);
            maxY = Math.max(maxY, pt[1]);
        }
        return (i: number) => {
            const x = data1[i];
            const y = data2[i];
            if (x < minX || x > maxX || y < minY || y > maxY || !Number.isFinite(x) || !Number.isFinite(y)) {
                return false;
            }
            if (isRegionScoped) {
                if (!regionData || wantedRegionIndex === -1) return false;
                if (regionData[i] !== wantedRegionIndex) return false;
            }
            let inside = false;
            for (let curr = 0, prev = points.length - 1; curr < points.length; prev = curr++) {
                const xi = points[curr][0];
                const yi = points[curr][1];
                const xj = points[prev][0];
                const yj = points[prev][1];
                const intersect = yi > y !== yj > y && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi;
                if (intersect) inside = !inside;
            }
            return inside;
        };
    },
};

const N = 4000;
const SPECIAL = [Number.NaN, Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY, 0, -0, -3, 3, 0.1];

/** Values on a coarse grid (so bounds and vertices land exactly on data) plus specials. */
function randomColumn(rand: () => number, size = N) {
    const out = new Float32Array(size);
    for (let i = 0; i < size; i++) {
        const r = rand();
        out[i] = r < 0.08 ? SPECIAL[Math.floor(rand() * SPECIAL.length)] : Math.round(rand() * 40 - 20) / 2;
    }
    return out;
}
const randomBound = (rand: () => number) =>
    rand() < 0.1 ? SPECIAL[Math.floor(rand() * SPECIAL.length)] : Math.round(rand() * 40 - 20) / 2;

function randomPolygon(rand: () => number) {
    const n = 3 + Math.floor(rand() * 7);
    // sometimes entirely negative, which exercises the Number.MIN_VALUE bounding-box seeds
    const offset = rand() < 0.3 ? -15 : 0;
    return Array.from({ length: n }, () => [randomBound(rand) + offset, randomBound(rand) + offset]) as [
        number,
        number,
    ][];
}

function mask(spec: FilterSpec, size = N) {
    const out = new Uint8Array(size);
    evaluateFilterSpec(spec, size, out);
    return out;
}
/** Compares a mask against a predicate, reporting the first row that differs. */
function expectMatches(spec: FilterSpec, keep: (i: number) => boolean) {
    const m = mask(spec);
    const row = m.findIndex((v, i) => v !== (keep(i) ? 1 : 0));
    expect(row, `first differing row (mask ${m[row]})`).toBe(-1);
}

describe("evaluateFilterSpec matches the legacy predicates", () => {
    test.each([1, 2, 3, 4, 5, 6])("range, seed %i", (seed) => {
        const rand = rng(seed);
        for (const data of [randomColumn(rand), Int32Array.from(randomColumn(rand), (v) => v | 0)]) {
            for (let k = 0; k < 40; k++) {
                const [min, max] = [randomBound(rand), randomBound(rand)];
                expectMatches({ kind: "range", data, min, max }, legacy.range(data, min, max));
            }
        }
    });

    test.each([1, 2, 3, 4, 5, 6])("square, seed %i", (seed) => {
        const rand = rng(seed);
        const [data1, data2] = [randomColumn(rand), randomColumn(rand)];
        for (let k = 0; k < 40; k++) {
            const range1: [number, number] = [randomBound(rand), randomBound(rand)];
            const range2: [number, number] = [randomBound(rand), randomBound(rand)];
            expectMatches(
                { kind: "square", data1, data2, range1, range2 },
                legacy.square(data1, data2, range1, range2),
            );
        }
    });

    test.each([1, 2, 3, 4, 5, 6])("poly, seed %i", (seed) => {
        const rand = rng(seed);
        const [data1, data2] = [randomColumn(rand), randomColumn(rand)];
        const regionData = Uint8Array.from({ length: N }, () => Math.floor(rand() * 3));
        for (let k = 0; k < 40; k++) {
            const points = randomPolygon(rand);
            expectMatches({ kind: "poly", data1, data2, points }, legacy.poly(data1, data2, points, false, null, -1));
            // region found, value missing from the column, column missing
            for (const [data, index] of [
                [regionData, Math.floor(rand() * 3)],
                [regionData, -1],
                [null, 0],
            ] as const) {
                expectMatches(
                    { kind: "poly", data1, data2, points, region: { data, index } },
                    legacy.poly(data1, data2, points, true, data, index),
                );
            }
        }
    });

    test("only rows below size are written", () => {
        const out = new Uint8Array(6).fill(7);
        evaluateFilterSpec({ kind: "range", data: new Float32Array([1, 2, 3, 4, 5, 6]), min: 0, max: 9 }, 4, out);
        expect(Array.from(out)).toEqual([1, 1, 1, 1, 7, 7]);
    });

    test("a short output is rejected", () => {
        expect(() => evaluateFilterSpec({ kind: "range", data: [1, 2], min: 0, max: 1 }, 2, new Uint8Array(1))).toThrow(
            /1 rows, expected 2/,
        );
    });
});

let ds: DataStore;
let workers: InProcessFilterWorker[];
let events: string[];
const restoreMinRows = RangeDimension.asyncFilterMinRows;

function store(seed: number, size = 2000) {
    const rand = rng(seed);
    const s = createStore(size, [
        { field: "x", datatype: "double", data: Array.from(randomColumn(rand, size)) },
        { field: "y", datatype: "double", data: Array.from(randomColumn(rand, size)) },
        { field: "r", datatype: "text", data: Array.from({ length: size }, () => (rand() < 0.5 ? "r1" : "r2")) },
        { field: "d", datatype: "text", data: Array.from({ length: size }, () => (rand() < 0.5 ? "p" : "q")) },
    ]);
    s.addListener("spy", (type: string) => events.push(type));
    return s;
}

beforeEach(() => {
    workers = stubWorker({ rangeFilter: true }).rangeFilterWorkers;
    vi.spyOn(console, "log").mockImplementation(() => {});
    events = [];
});
afterEach(async () => {
    if (ds) {
        await settleFilters(ds);
        assertConsistent(ds);
    }
    RangeDimension.asyncFilterMinRows = restoreMinRows;
    vi.restoreAllMocks();
    vi.unstubAllGlobals();
});

describe("RangeDimension async path", () => {
    const filters: [string, string[], object][] = [
        ["filterRange", ["x"], { min: -3, max: 4.5 }],
        ["filterSquare", ["x", "y"], { range1: [-5, 5], range2: [0, 10] }],
        [
            "filterPoly",
            ["x", "y"],
            [
                [-8, -8],
                [9, -2],
                [3, 9],
                [-2, 1],
            ],
        ],
        [
            "filterPoly",
            ["x", "y", "r"],
            {
                points: [
                    [-8, -8],
                    [9, -2],
                    [3, 9],
                ],
                regionField: "r",
                regionValue: "r2",
            },
        ],
    ];

    /** Applies the filter (under a background filter) and returns the committed local states. */
    async function run(method: string, columns: string[], args: object, forceAsync: boolean) {
        RangeDimension.asyncFilterMinRows = forceAsync ? 0 : Number.POSITIVE_INFINITY;
        ds = store(7);
        const dim = ds.getDimension("range_dimension");
        dim.setBackgroundFilter("d", "p");
        await dim.filter(method, columns, args);
        return { states: Array.from(dim.filterArray), passing: passingRows(ds) };
    }

    test.each(filters)("%s commits the same states as the sync path", async (method, columns, args) => {
        const sync = await run(method, columns, args, false);
        expect(workers).toHaveLength(0);
        const async = await run(method, columns, args, true);
        expect(workers[0].requests).toHaveLength(1);
        expect(async).toEqual(sync);
        expect(sync.passing.length).toBeGreaterThan(0);
    });

    test("filterSquare with a typed array for the first column", async () => {
        const xs = new Float32Array(2000).map((_, i) => i % 10);
        const sync = await run("filterSquare", [xs as never, "y"], { range1: [2, 4], range2: [-5, 5] }, false);
        const async = await run("filterSquare", [xs as never, "y"], { range1: [2, 4], range2: [-5, 5] }, true);
        expect(async).toEqual(sync);
    });

    test("below the threshold, filters commit synchronously without a worker", () => {
        RangeDimension.asyncFilterMinRows = 2001;
        ds = store(1);
        const dim = ds.getDimension("range_dimension");
        void dim.filter("filterRange", ["x"], { min: 0, max: 5 });
        expect(dim.pendingFilter).toBeNull();
        expect(ds.filterSize).toBeLessThan(2000);
        expect(workers).toHaveLength(0);
    });

    test("at the threshold, the result lands after filter() returns", async () => {
        RangeDimension.asyncFilterMinRows = 2000;
        ds = store(1);
        const dim = ds.getDimension("range_dimension");
        const done = dim.filter("filterRange", ["x"], { min: 0, max: 5 });
        expect(dim.pendingFilter).not.toBeNull();
        expect(ds.filterSize).toBe(2000);
        await expect(done).resolves.toBe(true);
        expect(ds.filterSize).toBeLessThan(2000);
        expect(events).toEqual(["filtered"]);
    });
});

describe("rapid brushing", () => {
    beforeEach(() => {
        RangeDimension.asyncFilterMinRows = 0;
        ds = store(3);
    });

    test("only the first and latest of a burst reach the worker, and only the latest commits", async () => {
        const dim = ds.getDimension("range_dimension");
        const results = Array.from({ length: 10 }, (_, k) => dim.filter("filterRange", ["x"], { min: -k, max: k }));
        expect(workers[0].requests).toHaveLength(1);

        expect(await Promise.all(results)).toEqual([...Array(9).fill(false), true]);
        expect(workers[0].requests.map((r) => r.spec)).toMatchObject([
            { min: -0, max: 0 },
            { min: -9, max: 9 },
        ]);
        expect(events).toEqual(["filtered"]);

        const x = ds.columnIndex.x?.data;
        if (!x) throw new Error("x not loaded");
        expect(passingRows(ds)).toEqual(
            Array.from({ length: 2000 }, (_, i) => i).filter((i) => x[i] >= -9 && x[i] <= 9),
        );
    });

    test("a brush that lands between results is still applied", async () => {
        const dim = ds.getDimension("range_dimension");
        await dim.filter("filterRange", ["x"], { min: 0, max: 1 });
        const second = dim.filter("filterSquare", ["x", "y"], { range1: [0, 3], range2: [0, 3] });
        await expect(second).resolves.toBe(true);
        expect(workers[0].requests).toHaveLength(2);
        expect(events).toEqual(["filtered", "filtered"]);
    });

    test("removeFilter drops the queued request before it reaches the worker", async () => {
        const dim = ds.getDimension("range_dimension");
        const first = dim.filter("filterRange", ["x"], { min: 0, max: 1 });
        const second = dim.filter("filterRange", ["x"], { min: 0, max: 2 });
        dim.removeFilter();
        expect(await Promise.all([first, second])).toEqual([false, false]);
        expect(workers[0].requests).toHaveLength(1);
        expect(ds.filterSize).toBe(2000);
    });

    test("destroy settles outstanding requests and stops the worker", async () => {
        const dim = ds.getDimension("range_dimension");
        const first = dim.filter("filterRange", ["x"], { min: 0, max: 1 });
        const second = dim.filter("filterRange", ["x"], { min: 0, max: 2 });
        dim.destroy();
        expect(await Promise.all([first, second])).toEqual([false, false]);
        expect(workers[0].terminated).toBe(true);
        expect(ds.filterSize).toBe(2000);
    });

    test("observable (proxied) polygon points are copied before posting", async () => {
        const dim = ds.getDimension("range_dimension");
        const square = [
            [-2, -2],
            [2, -2],
            [2, 2],
            [-2, 2],
        ].map((pt) => new Proxy(pt, {}));
        await expect(dim.filter("filterPoly", ["x", "y"], new Proxy(square, {}))).resolves.toBe(true);
        expect(workers[0].requests).toHaveLength(1);
        expect(ds.filterSize).toBeGreaterThan(0);
        expect(ds.filterSize).toBeLessThan(2000);
    });

    test("a spec that can't be posted is evaluated on the main thread", async () => {
        vi.spyOn(console, "warn").mockImplementation(() => {});
        const worker = new InProcessFilterWorker();
        const client = new RangeFilterClient(() => worker as never);
        const data = new Proxy(new Float32Array([1, 5, 3, Number.NaN]), {});
        const mask = await client.evaluate({ kind: "range", data, min: 2, max: 5 }, 4);
        expect(Array.from(mask ?? [])).toEqual([0, 1, 1, 0]);
        expect(worker.requests).toHaveLength(0);
    });

    test("if the worker fails, filtering falls back to the main thread", async () => {
        vi.spyOn(console, "warn").mockImplementation(() => {});
        const dim = ds.getDimension("range_dimension");
        const pending = dim.filter("filterRange", ["x"], { min: 0, max: 1 });
        workers[0].terminated = true;
        workers[0].onerror?.({ message: "boom", preventDefault() {} });
        await expect(pending).resolves.toBe(true);
        const afterFailure = ds.filterSize;
        expect(afterFailure).toBeLessThan(2000);

        await expect(dim.filter("filterRange", ["x"], { min: 0, max: 5 })).resolves.toBe(true);
        expect(ds.filterSize).toBeGreaterThan(afterFailure);
        expect(workers[0].requests).toHaveLength(1);
    });
});
