import { observable } from "mobx";
import { describe, expect, test } from "vitest";
import { compileCategoryFilter } from "@/lib/categoryFilterUtils";
import {
    getRowMultitextItems,
    getRowMultitextValueIndices,
    inferMultitextEmptyItem,
    MULTITEXT_EMPTY_INDEX,
} from "@/lib/multitext";
import { isArray } from "@/lib/utils";
import { rng } from "../datastore/dimensionHarness";

type Column = Parameters<typeof compileCategoryFilter>[0];

/** The per-row matcher compileCategoryFilter replaced, verbatim. */
function legacyRowMatches(rowIndex: number, column: Column, category: string | string[] | null) {
    if (!category || !column.data) {
        return false;
    }
    const selectedValues = (isArray(category) ? category : [category]).filter(
        (value): value is string => typeof value === "string" && value.length > 0,
    );
    if (selectedValues.length === 0) {
        return false;
    }
    if (column.datatype === "multitext") {
        const emptyItem = inferMultitextEmptyItem(column);
        const selectedValueIndices = new Set(
            selectedValues.map((value) => column.values.indexOf(value)).filter((valueIndex) => valueIndex !== -1),
        );
        const rowValueIndices = getRowMultitextValueIndices(column, rowIndex);
        if (rowValueIndices.some((valueIndex) => selectedValueIndices.has(valueIndex))) {
            return true;
        }
        const rowItems = getRowMultitextItems(column, rowIndex, { emptyItem });
        return selectedValues.some((value) => rowItems.includes(value));
    }
    const categoryValueIndices = selectedValues
        .map((value) => column.values.indexOf(value))
        .filter((valueIndex) => valueIndex !== -1);
    if (categoryValueIndices.length === 0) {
        return false;
    }
    return categoryValueIndices.includes(column.data[rowIndex]);
}

const ROWS = 300;
const pick = <T>(rand: () => number, xs: readonly T[]) => xs[Math.floor(rand() * xs.length)];

function columns(rand: () => number): [string, Column][] {
    // "a" twice: only its first index is the selected raw value
    const textValues = ["a", "b", "c", "a", "d e"];
    const listValues = ["a, b", "b", "c;d", "N/A", "", "b , e"];
    const packedValues = ["N/A", "a", "b", "c", "d"];
    const packed = new Uint16Array(ROWS * 3).fill(MULTITEXT_EMPTY_INDEX);
    for (let row = 0; row < ROWS; row++) {
        const n = Math.floor(rand() * 4);
        for (let k = 0; k < n; k++) packed[row * 3 + k] = Math.floor(rand() * packedValues.length);
    }
    return [
        [
            "text",
            {
                datatype: "text",
                values: textValues,
                data: Uint8Array.from({ length: ROWS }, () => Math.floor(rand() * (textValues.length + 1))),
            },
        ],
        [
            "multitext",
            {
                datatype: "multitext",
                values: listValues,
                stringLength: 1,
                data: Uint8Array.from({ length: ROWS }, () => Math.floor(rand() * listValues.length)),
            },
        ],
        [
            "packed multitext",
            { datatype: "multitext", values: packedValues, stringLength: 3, delimiter: ";", data: packed },
        ],
    ];
}

const CANDIDATES = ["a", "b", "c", "d", "e", "d e", "a, b", "c;d", "N/A", "", "zzz"];

describe("compileCategoryFilter matches the per-row matcher it replaced", () => {
    test.each([1, 2, 3, 4])("seed %i", (seed) => {
        const rand = rng(seed);
        for (const [name, column] of columns(rand)) {
            for (let k = 0; k < 30; k++) {
                const roll = rand();
                const category: string | string[] | null =
                    roll < 0.1 ? null : roll < 0.4 ? pick(rand, CANDIDATES) : CANDIDATES.filter(() => rand() < 0.3);
                const matches = compileCategoryFilter(column, category);
                for (let row = 0; row < ROWS; row++) {
                    expect(matches(row), `${name} ${JSON.stringify(category)} row ${row}`).toBe(
                        legacyRowMatches(row, column, category),
                    );
                }
            }
        }
    });

    test("an observable category is read once, when compiling", () => {
        const category = observable(["b"]);
        const [, column] = columns(rng(1))[0];
        const matches = compileCategoryFilter(column, category);
        const before = Array.from({ length: ROWS }, (_, row) => matches(row));
        category.push("c");
        expect(Array.from({ length: ROWS }, (_, row) => matches(row))).toEqual(before);
        expect(before.some(Boolean)).toBe(true);
    });
});
