import { describe, expect, test } from "vitest";
import { ownSelectionMask } from "@/react/scatter_state";

describe("telling a chart's own selection apart", () => {
    const rows = new Uint32Array([0, 1, 2, 3, 4]);
    const owner = (local: number[], filterMethod?: string) => ({
        getLocalFilter: () => new Uint8Array(local),
        filterMethod,
    });

    test("there is nothing to fade without a filter of its own", () => {
        expect(ownSelectionMask(rows, 5, owner([1, 1, 0, 0, 0]))).toBeNull();
    });

    test("marks the drawn rows its own filter leaves out", () => {
        // row 5 is not drawn (another chart's filter), so it is not marked either
        const mask = ownSelectionMask(rows, 6, owner([0, 1, 0, 3, 2, 1], "filterScreenPoly"));
        expect(Array.from(mask ?? [])).toEqual([0, 1, 0, 1, 0, 0]);
    });

    test("there is nothing to fade when the filter leaves no drawn row out", () => {
        expect(ownSelectionMask(rows, 5, owner([0, 0, 0, 0, 2], "filterScreenPoly"))).toBeNull();
    });
});
