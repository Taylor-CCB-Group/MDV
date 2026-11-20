import { describe, expect, test } from "vitest";
import { splitBySelection } from "@/react/scatter_state";

describe("telling a chart's own selection apart", () => {
    const rows = new Uint32Array([0, 1, 2, 3, 4]);
    const owner = (local: number[], filterMethod?: string) => ({
        getLocalFilter: () => new Uint8Array(local),
        filterMethod,
    });

    test("every row is selected without a filter of its own", () => {
        const split = splitBySelection(rows, owner([1, 1, 0, 0, 0]));
        expect(split.selectedRows).toBe(rows);
        expect(split.unselectedRows).toHaveLength(0);
    });

    test("rows its own filter leaves out are the unselected ones", () => {
        const split = splitBySelection(rows, owner([0, 1, 0, 3, 2], "filterPoly"));
        expect([...split.selectedRows]).toEqual([0, 2, 4]);
        expect([...split.unselectedRows]).toEqual([1, 3]);
    });

    test("keeps the rows as they were when the filter leaves none out", () => {
        const split = splitBySelection(rows, owner([0, 0, 0, 0, 2], "filterPoly"));
        expect(split.selectedRows).toBe(rows);
    });
});
