import Dimension from "@/datastore/Dimension";
import { describe, expect, test } from "vitest";

function createParent(size: number, bgData?: Uint8Array) {
    return {
        size,
        filterArray: new Uint8Array(size),
        filterSize: size,
        columnIndex: {
            group: {
                values: ["A", "B"],
                data: bgData ?? new Uint8Array(size),
            },
        },
        _callListeners: () => {},
        dimensions: [] as unknown[],
    };
}

function createDimension(size: number, bgData?: Uint8Array) {
    const parent = createParent(size, bgData);
    const dim = new Dimension(parent as never);
    parent.dimensions.push(dim);
    return { parent, dim };
}

function assertInvariant(parent: ReturnType<typeof createParent>, dim: Dimension) {
    let unfiltered = 0;
    for (let i = 0; i < parent.size; i++) {
        const state = dim.filterArray[i];
        const isFiltered = state === 1 || state === 3;
        expect(parent.filterArray[i]).toBe(isFiltered ? 1 : 0);
        if (!isFiltered) {
            unfiltered++;
        }
    }
    expect(parent.filterSize).toBe(unfiltered);
}

describe("Dimension 0/1/2/3 background-filter transitions", () => {
    test("brush only mutates visible states and keeps global counts consistent", () => {
        const { parent, dim } = createDimension(4);
        dim.filterArray.set([0, 1, 2, 3]);
        parent.filterArray.set([0, 1, 0, 1]);
        parent.filterSize = 2;

        const includeRows = new Set([0, 1]);
        dim.filterPredicate({ predicate: (i: number) => includeRows.has(i) }, []);
        // background-hidden rows promoted to state 3 when local filter active
        expect(Array.from(dim.filterArray)).toEqual([0, 0, 3, 3]);
        assertInvariant(parent, dim);

        const includeOnlyFirst = new Set([0]);
        dim.filterPredicate(
            { predicate: (i: number) => includeOnlyFirst.has(i) },
            [],
        );

        // exclude row 1: 0->1
        expect(Array.from(dim.filterArray)).toEqual([0, 1, 3, 3]);
        assertInvariant(parent, dim);
    });

    test("set/clear background transitions avoid count drift", () => {
        const bg = new Uint8Array([0, 1, 0, 1]); // A, B, A, B
        const { parent, dim } = createDimension(4, bg);
        dim.filterArray.set([0, 1, 0, 0]);
        parent.filterArray.set([0, 1, 0, 0]);
        parent.filterSize = 3;

        dim.setBackgroundFilter("group", ["B"]);
        expect(Array.from(dim.filterArray)).toEqual([2, 1, 2, 0]);
        assertInvariant(parent, dim);

        dim.clearBackgroundFilter();
        expect(Array.from(dim.filterArray)).toEqual([0, 1, 0, 0]);
        assertInvariant(parent, dim);
    });

    test("removeFilter clears only local contribution", () => {
        const { parent, dim } = createDimension(4);
        dim.filterArray.set([0, 1, 2, 3]);
        parent.filterArray.set([0, 1, 0, 1]);
        parent.filterSize = 2;
        dim.filterMethod = "filterPredicate";

        dim.removeFilter(false);

        expect(Array.from(dim.filterArray)).toEqual([0, 0, 2, 2]);
        assertInvariant(parent, dim);
    });
});
