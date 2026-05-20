import CategoryDimension from "@/datastore/CategoryDimension";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";

function createParent(column: any, size: number) {
    return {
        size,
        filterArray: new Uint8Array(size),
        filterSize: size,
        columnIndex: {
            tags: column,
        },
        _callListeners: vi.fn(),
    };
}

function getIncludedRows(dimension: CategoryDimension, size: number) {
    const rows: number[] = [];
    for (let index = 0; index < size; index += 1) {
        if (dimension.filterArray[index] === 0) {
            rows.push(index);
        }
    }
    return rows;
}

describe("CategoryDimension multitext filtering", () => {
    beforeEach(() => {
        const WorkerMock = vi
            .fn()
            .mockImplementation(function MockWorker() {
                return {
                    addEventListener: vi.fn(),
                    postMessage: vi.fn(),
                    terminate: vi.fn(),
                };
            });
        vi.stubGlobal(
            "Worker",
            WorkerMock,
        );
    });

    afterEach(() => {
        vi.unstubAllGlobals();
    });

    test("matches legacy combined values by individual items", () => {
        const column = {
            field: "tags",
            name: "tags",
            datatype: "multitext" as const,
            delimiter: ",",
            stringLength: 1,
            values: ["a", "b", "a, b"],
            data: new Uint16Array([0, 1, 2]),
        };
        const parent = createParent(column, 3);
        const dimension = new CategoryDimension(parent as never);
        const filter = ["b"] as string[] & { operand?: "or" | "and" };
        filter.operand = "or";

        dimension.filterCategories(filter, ["tags"]);

        expect(getIncludedRows(dimension, 3)).toEqual([1, 2]);
    });

    test("supports AND matching across combined legacy values", () => {
        const column = {
            field: "tags",
            name: "tags",
            datatype: "multitext" as const,
            delimiter: ",",
            stringLength: 1,
            values: ["a", "b", "a, b", "b, c"],
            data: new Uint16Array([0, 1, 2, 3]),
        };
        const parent = createParent(column, 4);
        const dimension = new CategoryDimension(parent as never);
        const filter = ["a", "b"] as string[] & { operand?: "or" | "and" };
        filter.operand = "and";

        dimension.filterCategories(filter, ["tags"]);

        expect(getIncludedRows(dimension, 4)).toEqual([2]);
    });

    test("matches packed multitext rows by individual items", () => {
        const column = {
            field: "tags",
            name: "tags",
            datatype: "multitext" as const,
            delimiter: ",",
            stringLength: 3,
            values: ["N/A", "a", "b"],
            data: new Uint16Array([
                0, 65535, 65535,
                2, 65535, 65535,
                1, 2, 65535,
            ]),
        };
        const parent = createParent(column, 3);
        const dimension = new CategoryDimension(parent as never);
        const filter = ["b"] as string[] & { operand?: "or" | "and" };
        filter.operand = "or";

        dimension.filterCategories(filter, ["tags"]);

        expect(getIncludedRows(dimension, 3)).toEqual([1, 2]);
    });

    test("can require an exact single-item multitext match", () => {
        const column = {
            field: "tags",
            name: "tags",
            datatype: "multitext" as const,
            delimiter: ",",
            stringLength: 1,
            values: ["a", "b", "a, b"],
            data: new Uint16Array([0, 1, 2]),
        };
        const parent = createParent(column, 3);
        const dimension = new CategoryDimension(parent as never);

        // Use a string filter instead of an object
        dimension.filterCategories("a", ["tags"]);

        expect(getIncludedRows(dimension, 3)).toEqual([0]);
    });
});
