import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import type DataStore from "@/datastore/DataStore";
import type Dimension from "@/datastore/Dimension";
import CategoryDimension from "@/datastore/CategoryDimension";
import { assertConsistent, createStore, passingRows, rng, stubWorker } from "./dimensionHarness";

/**
 * Characterization tests for the Dimension filtering engine, run against a real DataStore.
 * Every test finishes with `assertConsistent`, which checks the parent ref-counts against the
 * per-dimension 0/1/2/3 states.
 */

// row:          0    1    2    3    4    5    6    7
const X = [0, 1, 2, 3, 4, 5, Number.NaN, 7];
const Y = [0, 1, 4, 9, 16, 25, 36, Number.NaN];
const C = ["a", "b", "a", "c", "b", "a", "c", "b"];
const D = ["p", "p", "q", "q", "p", "q", "p", "q"];
const R = ["r1", "r1", "r1", "r1", "r2", "r2", "r2", "r2"];

function fixture() {
    return createStore(8, [
        { field: "x", datatype: "double", data: X },
        { field: "y", datatype: "double", data: Y },
        { field: "c", datatype: "text", data: C },
        { field: "d", datatype: "text", data: D },
        { field: "r", datatype: "text", data: R },
    ]);
}

const ALL = [0, 1, 2, 3, 4, 5, 6, 7];

let ds: DataStore;
beforeEach(() => {
    stubWorker();
    ds = fixture();
});
afterEach(() => {
    assertConsistent(ds);
    vi.unstubAllGlobals();
});

describe("RangeDimension", () => {
    test("filterRange keeps inclusive bounds and drops NaN", () => {
        const dim = ds.getDimension("range_dimension");
        dim.filter("filterRange", ["x"], { min: 1, max: 4 });
        expect(passingRows(ds)).toEqual([1, 2, 3, 4]);
        expect(ds.filterSize).toBe(4);
    });

    test("re-filtering narrows and widens without count drift", () => {
        const dim = ds.getDimension("range_dimension");
        dim.filter("filterRange", ["x"], { min: 1, max: 4 });
        dim.filter("filterRange", ["x"], { min: 3, max: 7 });
        expect(passingRows(ds)).toEqual([3, 4, 5, 7]);
        assertConsistent(ds);
        dim.filter("filterRange", ["x"], { min: -100, max: 100 });
        expect(passingRows(ds)).toEqual([0, 1, 2, 3, 4, 5, 7]);
    });

    test("filterSquare requires both axes in range and finite", () => {
        const dim = ds.getDimension("range_dimension");
        dim.filter("filterSquare", ["x", "y"], { range1: [1, 7], range2: [1, 16] });
        expect(passingRows(ds)).toEqual([1, 2, 3, 4]);
    });

    test("filterSquare accepts a typed array in place of the first column", () => {
        const dim = ds.getDimension("range_dimension");
        const xs = new Float32Array([5, 5, 5, 5, 0, 0, 0, 0]);
        dim.filter("filterSquare", [xs, "y"] as never, { range1: [4, 6], range2: [-1, 100] });
        expect(passingRows(ds)).toEqual([0, 1, 2, 3]);
    });

    test("filterPoly keeps rows inside the polygon", () => {
        const dim = ds.getDimension("range_dimension");
        const square = [
            [0.5, 0.5],
            [4.5, 0.5],
            [4.5, 10],
            [0.5, 10],
        ];
        dim.filter("filterPoly", ["x", "y"], square);
        expect(passingRows(ds)).toEqual([1, 2, 3]);
    });

    test("region-scoped filterPoly only keeps rows in that region", () => {
        const dim = ds.getDimension("range_dimension");
        const points = [
            [0.5, 0.5],
            [5.5, 0.5],
            [5.5, 30],
            [0.5, 30],
        ];
        dim.filter("filterPoly", ["x", "y", "r"], { points, regionField: "r", regionValue: "r1" });
        expect(passingRows(ds)).toEqual([1, 2, 3]);
        dim.filter("filterPoly", ["x", "y", "r"], { points, regionField: "r", regionValue: "r2" });
        expect(passingRows(ds)).toEqual([4, 5]);
        dim.filter("filterPoly", ["x", "y", "r"], { points, regionField: "r", regionValue: "nope" });
        expect(passingRows(ds)).toEqual([]);
    });
});

describe("CategoryDimension", () => {
    test("filterCategories with a single category string", () => {
        const dim = ds.getDimension("category_dimension");
        dim.filter("filterCategories", ["c"], "a");
        expect(passingRows(ds)).toEqual([0, 2, 5]);
    });

    test("filterCategories with a list of categories", () => {
        const dim = ds.getDimension("category_dimension");
        dim.filter("filterCategories", ["c"], ["a", "c"]);
        expect(passingRows(ds)).toEqual([0, 2, 3, 5, 6]);
        dim.filter("filterCategories", ["c"], ["b"]);
        expect(passingRows(ds)).toEqual([1, 4, 7]);
    });

    test("filterCategories with an unknown category excludes everything", () => {
        const dim = ds.getDimension("category_dimension");
        dim.filter("filterCategories", ["c"], ["zz"]);
        expect(passingRows(ds)).toEqual([]);
    });

    test.each([
        ["a", [0, 2, 5]],
        [["a"], [0, 2, 5]],
        [["a", "c"], [0, 2, 3, 5, 6]],
        [["a", "b", "c"], ALL],
    ])("filterCategoriesPredicate(%j) matches filterCategories", (cats, expected) => {
        const dim = ds.getDimension("category_dimension");
        dim.filter("filterCategoriesPredicate", ["c"], cats);
        expect(passingRows(ds)).toEqual(expected);
    });

    test("filterMultipleCategories requires every column to match", () => {
        const dim = ds.getDimension("category_dimension");
        dim.filter("filterMultipleCategories", ["c", "d"], ["a", "q"]);
        expect(passingRows(ds)).toEqual([2, 5]);
    });

    test("filterCategoryPair requires both columns to match", () => {
        const dim = ds.getDimension("category_dimension");
        dim.filter("filterCategoryPair", ["c", "d"], ["a", "q"]);
        expect(passingRows(ds)).toEqual([2, 5]);
        dim.filter("filterCategoryPair", ["c", "d"], ["a", "nope"]);
        expect(passingRows(ds)).toEqual([]);
    });

    test("filterUnique matches every space-separated search term", () => {
        const dim = ds.getDimension("category_dimension");
        dim.filter("filterUnique", ["c"], "a");
        expect(passingRows(ds)).toEqual([0, 2, 5]);
    });
});

describe("ValueSetDimension", () => {
    test("filters text columns by value", () => {
        const dim = ds.getDimension("valueset_dimension");
        dim.filter("filterValueset", ["c"], new Set(["b", "c"]));
        expect(passingRows(ds)).toEqual([1, 3, 4, 6, 7]);
    });

    test("filters numeric columns by raw value", () => {
        const dim = ds.getDimension("valueset_dimension");
        dim.filter("filterValueset", ["x"], new Set([2, 4]));
        expect(passingRows(ds)).toEqual([2, 4]);
    });

    test("an empty set leaves the previous filter in place", () => {
        const dim = ds.getDimension("valueset_dimension");
        dim.filter("filterValueset", ["c"], new Set(["b"]));
        dim.filter("filterValueset", ["c"], new Set());
        expect(passingRows(ds)).toEqual([1, 4, 7]);
    });
});

describe("CatColDimension", () => {
    test("filterCatCol keeps the category where the value exceeds the threshold", () => {
        const dim = ds.getDimension("catcol_dimension");
        dim.filter("filterCatCol", ["c", "x"], { cat: "a", threshold: 1 });
        expect(passingRows(ds)).toEqual([2, 5]);
    });
});

describe("base Dimension", () => {
    test("filterPredicate via filter()", () => {
        const dim = ds.getDimension("base_dimension");
        dim.filter("filterPredicate", [], { predicate: (i: number) => i % 3 === 0 });
        expect(passingRows(ds)).toEqual([0, 3, 6]);
    });

    test("filterOnIndex keeps only the given indices", () => {
        const dim = ds.getDimension("base_dimension");
        dim.filter("filterOnIndex", [], new Set([1, 3]));
        expect(passingRows(ds)).toEqual([1, 3]);
        dim.filter("filterOnIndex", [], new Set([3, 4]));
        expect(passingRows(ds)).toEqual([3, 4]);
    });
});

describe("multiple dimensions", () => {
    test("passing rows are the intersection; removing one restores only its exclusions", () => {
        const range = ds.getDimension("range_dimension");
        const cat = ds.getDimension("category_dimension");
        range.filter("filterRange", ["x"], { min: 1, max: 4 });
        cat.filter("filterCategories", ["c"], ["a"]);
        expect(passingRows(ds)).toEqual([2]);
        assertConsistent(ds);
        cat.removeFilter();
        expect(passingRows(ds)).toEqual([1, 2, 3, 4]);
        assertConsistent(ds);
        range.removeFilter();
        expect(passingRows(ds)).toEqual(ALL);
    });

    test("removeAllFilters clears everything except noClear dimensions", () => {
        const range = ds.getDimension("range_dimension");
        const cat = ds.getDimension("category_dimension");
        const sticky = ds.getDimension("valueset_dimension");
        sticky.setNoClear(true);
        range.filter("filterRange", ["x"], { min: 1, max: 4 });
        cat.filter("filterCategories", ["c"], ["a"]);
        sticky.filter("filterValueset", ["d"], new Set(["q"]));
        expect(passingRows(ds)).toEqual([2]);

        ds.removeAllFilters();

        expect(passingRows(ds)).toEqual([2, 3, 5, 7]);
        expect(range.filterMethod).toBeNull();
        expect(cat.filterMethod).toBeNull();
        expect(sticky.filterMethod).toBe("filterValueset");
    });

    test("a dimension can filter again after removeAllFilters", () => {
        const range = ds.getDimension("range_dimension");
        range.filter("filterRange", ["x"], { min: 1, max: 4 });
        ds.removeAllFilters();
        range.filter("filterRange", ["x"], { min: 3, max: 5 });
        expect(passingRows(ds)).toEqual([3, 4, 5]);
    });

    test("destroy removes the dimension's contribution and is idempotent", () => {
        const range = ds.getDimension("range_dimension");
        const cat = ds.getDimension("category_dimension");
        range.filter("filterRange", ["x"], { min: 1, max: 4 });
        cat.filter("filterCategories", ["c"], ["a"]);
        cat.destroy();
        cat.destroy();
        expect(ds.dimensions).toEqual([range]);
        expect(passingRows(ds)).toEqual([1, 2, 3, 4]);
    });
});

describe("notifications", () => {
    function listen() {
        const events: [string, unknown][] = [];
        ds.addListener("spy", (type: string, data: unknown) => events.push([type, data]));
        return events;
    }

    test("filter notifies with the dimension unless notify is false", () => {
        const events = listen();
        const dim = ds.getDimension("range_dimension");
        dim.filter("filterRange", ["x"], { min: 1, max: 4 });
        expect(events).toEqual([["filtered", dim]]);
        dim.filter("filterRange", ["x"], { min: 1, max: 2 }, false);
        expect(events).toHaveLength(1);
    });

    test("listeners observe the committed state", () => {
        const dim = ds.getDimension("range_dimension");
        let seen: number[] = [];
        ds.addListener("spy", () => {
            seen = passingRows(ds);
        });
        dim.filter("filterRange", ["x"], { min: 1, max: 4 });
        expect(seen).toEqual([1, 2, 3, 4]);
    });

    test("removeFilter notifies only when a filter was active", () => {
        const events = listen();
        const dim = ds.getDimension("range_dimension");
        dim.removeFilter();
        expect(events).toEqual([]);
        dim.filter("filterRange", ["x"], { min: 1, max: 4 }, false);
        dim.removeFilter();
        expect(events).toEqual([["filtered", dim]]);
    });

    test("removeAllFilters notifies once", () => {
        const events = listen();
        ds.getDimension("range_dimension").filter("filterRange", ["x"], { min: 1, max: 4 }, false);
        ds.getDimension("category_dimension").filter("filterCategories", ["c"], "a", false);
        ds.removeAllFilters();
        expect(events).toEqual([["filtered", "all_removed"]]);
    });
});

describe("data changes", () => {
    test("dataChanged re-evaluates only dimensions filtering on the changed columns", () => {
        const range = ds.getDimension("range_dimension");
        const cat = ds.getDimension("category_dimension");
        range.filter("filterRange", ["x"], { min: 1, max: 4 });
        cat.filter("filterCategories", ["c"], ["a", "b"]);
        const catMethod = vi.spyOn(CategoryDimension.prototype, "filterCategories");
        const events: [string, unknown][] = [];
        ds.addListener("spy", (type: string, data: unknown) => events.push([type, data]));

        const x = ds.columnIndex.x?.data;
        if (!x) throw new Error("x not loaded");
        x[0] = 2; // row 0 now in range
        ds.dataChanged(["x"], false);

        expect(passingRows(ds)).toEqual([0, 1, 2, 4]);
        expect(catMethod).not.toHaveBeenCalled();
        expect(events.map(([t]) => t)).toEqual(["filtered", "data_changed"]);
        expect(events[1][1]).toEqual({ columns: ["x"], hasFiltered: true });
    });

    test("dataChanged on unrelated columns does not refilter", () => {
        const range = ds.getDimension("range_dimension");
        range.filter("filterRange", ["x"], { min: 1, max: 4 });
        const events: string[] = [];
        ds.addListener("spy", (type: string) => events.push(type));
        ds.dataChanged(["c"], false);
        expect(events).toEqual(["data_changed"]);
    });

    test("addDataToStore extends local filters and re-evaluates over the new rows", () => {
        const range = ds.getDimension("range_dimension");
        range.filter("filterRange", ["x"], { min: 1, max: 4 });
        ds.addDataToStore(
            [
                { column: "x", data: new Float32Array([2, 9]).buffer },
                { column: "y", data: new Float32Array([0, 0]).buffer },
                { column: "c", data: new Uint8Array([0, 0]).buffer },
                { column: "d", data: new Uint8Array([0, 0]).buffer },
                { column: "r", data: new Uint8Array([0, 0]).buffer },
            ],
            2,
        );
        expect(ds.size).toBe(10);
        expect(range.filterArray.length).toBe(10);
        expect(passingRows(ds)).toEqual([1, 2, 3, 4, 8]);
    });
});

describe("background filters", () => {
    test("background-hidden rows stay globally visible until a local filter is applied", () => {
        const dim = ds.getDimension("range_dimension");
        dim.setBackgroundFilter("c", "a");
        expect(Array.from(dim.filterArray)).toEqual([0, 2, 0, 2, 2, 0, 2, 2]);
        expect(passingRows(ds)).toEqual(ALL);
    });

    test("a local filter also excludes background-hidden rows", () => {
        const dim = ds.getDimension("range_dimension");
        dim.setBackgroundFilter("c", "a");
        dim.filter("filterRange", ["x"], { min: 1, max: 4 });
        expect(Array.from(dim.filterArray)).toEqual([1, 3, 0, 3, 3, 1, 3, 3]);
        expect(passingRows(ds)).toEqual([2]);
        assertConsistent(ds);

        dim.removeFilter();
        expect(Array.from(dim.filterArray)).toEqual([0, 2, 0, 2, 2, 0, 2, 2]);
        expect(passingRows(ds)).toEqual(ALL);
    });

    test("setting a background filter while a local filter is active re-applies it", () => {
        const dim = ds.getDimension("range_dimension");
        dim.filter("filterRange", ["x"], { min: 1, max: 4 });
        dim.setBackgroundFilter("c", ["a", "b"]);
        expect(passingRows(ds)).toEqual([1, 2, 4]);
        assertConsistent(ds);
        dim.setBackgroundFilter("c", "c");
        expect(passingRows(ds)).toEqual([3]);
        assertConsistent(ds);
        dim.clearBackgroundFilter();
        expect(passingRows(ds)).toEqual([1, 2, 3, 4]);
    });

    test("an empty category list clears the background filter", () => {
        const dim = ds.getDimension("range_dimension");
        dim.setBackgroundFilter("c", "a");
        dim.setBackgroundFilter("c", []);
        expect(dim.bgfData).toBeNull();
        expect(Array.from(dim.filterArray)).toEqual(new Array(8).fill(0));
    });

    test("removeAllFilters keeps background states and drops local ones", () => {
        const dim = ds.getDimension("range_dimension");
        const other = ds.getDimension("category_dimension");
        dim.setBackgroundFilter("c", "a");
        dim.filter("filterRange", ["x"], { min: 1, max: 4 });
        other.filter("filterCategories", ["d"], "p");
        ds.removeAllFilters();
        expect(Array.from(dim.filterArray)).toEqual([0, 2, 0, 2, 2, 0, 2, 2]);
        expect(passingRows(ds)).toEqual(ALL);
        assertConsistent(ds);
        dim.filter("filterRange", ["x"], { min: 3, max: 5 });
        expect(passingRows(ds)).toEqual([5]);
    });

    test("background filters respect other dimensions' exclusions", () => {
        const dim = ds.getDimension("range_dimension");
        const other = ds.getDimension("category_dimension");
        other.filter("filterCategories", ["d"], "q");
        dim.setBackgroundFilter("c", "a");
        dim.filter("filterRange", ["x"], { min: 0, max: 7 });
        expect(passingRows(ds)).toEqual([2, 5]);
    });

    // Known defect: these methods write 0/1 directly instead of going through the 0/1/2/3 state
    // machine, so background-hidden rows are un-hidden (and the parent count drifts).
    test.fails.each([
        ["filterCategories", ["c"], ["a"]],
        ["filterMultipleCategories", ["c"], ["a"]],
    ])("%s honours an active background filter", (method, columns, args) => {
        const dim = ds.getDimension("category_dimension");
        dim.setBackgroundFilter("d", "p");
        dim.filter(method, columns, args);
        expect(passingRows(ds)).toEqual([0]);
        assertConsistent(ds);
    });

    test.fails("filterOnIndex honours an active background filter", () => {
        const dim = ds.getDimension("base_dimension");
        dim.setBackgroundFilter("d", "p");
        dim.filter("filterOnIndex", [], new Set([0, 1, 2]));
        expect(passingRows(ds)).toEqual([0, 1]);
        assertConsistent(ds);
    });
});

/**
 * Model-based test: applies a random sequence of operations to several dimensions and, after
 * each step, compares every local state against an independent model of the 0/1/2/3 rules.
 */
describe("randomized against a reference model", () => {
    type Model = { include: ((i: number) => boolean) | null; bg: Set<string> | null };
    type Op = (rand: () => number) => void;

    const pick = <T,>(rand: () => number, xs: readonly T[]) => xs[Math.floor(rand() * xs.length)];
    const subset = (rand: () => number, xs: readonly string[]) => {
        const s = xs.filter(() => rand() < 0.5);
        return s.length ? s : [xs[0]];
    };

    function expectedState(m: Model, i: number) {
        if (m.bg && !m.bg.has(D[i])) return m.include ? 3 : 2;
        return m.include && !m.include(i) ? 1 : 0;
    }

    test.each([1, 2, 3, 4, 5, 6, 7, 8])("seed %i", (seed) => {
        const rand = rng(seed);
        const dims: [Dimension, Model, Op][] = [];

        const range = ds.getDimension("range_dimension");
        const rangeModel: Model = { include: null, bg: null };
        dims.push([
            range,
            rangeModel,
            (r) => {
                const a = Math.floor(r() * 8) - 1;
                const b = a + Math.floor(r() * 6);
                range.filter("filterRange", ["x"], { min: a, max: b });
                rangeModel.include = (i) => X[i] >= a && X[i] <= b;
            },
        ]);

        const pred = ds.getDimension("category_dimension");
        const predModel: Model = { include: null, bg: null };
        dims.push([
            pred,
            predModel,
            (r) => {
                const cats = subset(r, ["a", "b", "c"]);
                pred.filter("filterCategoriesPredicate", ["c"], cats);
                predModel.include = (i) => cats.includes(C[i]);
            },
        ]);

        const vs = ds.getDimension("valueset_dimension");
        const vsModel: Model = { include: null, bg: null };
        dims.push([
            vs,
            vsModel,
            (r) => {
                const regions = subset(r, ["r1", "r2"]);
                vs.filter("filterValueset", ["r"], new Set(regions));
                vsModel.include = (i) => regions.includes(R[i]);
            },
        ]);

        for (let step = 0; step < 60; step++) {
            const roll = rand();
            if (roll < 0.04) {
                ds.removeAllFilters();
                for (const [, m] of dims) m.include = null;
            } else {
                const [dim, model, applyFilter] = pick(rand, dims);
                if (roll < 0.55) {
                    applyFilter(rand);
                } else if (roll < 0.7) {
                    dim.removeFilter();
                    model.include = null;
                } else if (roll < 0.88) {
                    const cats = subset(rand, ["p", "q"]);
                    dim.setBackgroundFilter("d", cats);
                    model.bg = new Set(cats);
                } else {
                    dim.clearBackgroundFilter();
                    model.bg = null;
                }
            }
            for (const [dim, model] of dims) {
                const expected = ALL.map((i) => expectedState(model, i));
                expect(Array.from(dim.filterArray), `step ${step}`).toEqual(expected);
            }
            assertConsistent(ds);
        }
    });
});
