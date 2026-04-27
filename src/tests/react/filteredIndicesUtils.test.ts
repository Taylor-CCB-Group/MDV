import { shouldRefreshFilteredIndices } from "@/react/filteredIndicesUtils";
import { describe, expect, test } from "vitest";

describe("shouldRefreshFilteredIndices", () => {
    test("refreshes for datastore events that invalidate filtered indices", () => {
        expect(shouldRefreshFilteredIndices("filtered")).toBe(true);
        expect(shouldRefreshFilteredIndices("data_changed")).toBe(true);
        expect(shouldRefreshFilteredIndices("data_added")).toBe(true);
    });

    test("ignores unrelated datastore events", () => {
        expect(shouldRefreshFilteredIndices("data_highlighted")).toBe(false);
        expect(shouldRefreshFilteredIndices("column_removed")).toBe(false);
    });
});
