import {
    getAvailableCategorySelectionValues,
    getCategorySelectionDropdownKey,
} from "@/react/components/categorySelectionUtils";
import { describe, expect, test, vi } from "vitest";

describe("categorySelectionUtils", () => {
    test("reads single-item values from multitext source columns", () => {
        const dataStore = {
            columnIndex: {
                __tags: {
                    datatype: "multitext",
                    delimiter: ",",
                    stringLength: 1,
                    values: ["a", "a, b", "b"],
                },
            },
            getColumnValues: vi.fn(),
        } as any;

        expect(getAvailableCategorySelectionValues(dataStore, "__tags")).toEqual([
            "a",
            "b",
        ]);
    });

    test("includes source column in dropdown remount key", () => {
        expect(getCategorySelectionDropdownKey("col_a", ["x", "y"])).not.toBe(
            getCategorySelectionDropdownKey("col_b", ["x", "y"]),
        );
    });
});
