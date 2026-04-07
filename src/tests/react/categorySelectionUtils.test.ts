import {
    getAvailableCategorySelectionValues,
    getCategorySelectionDropdownKey,
    shouldRefreshCategorySelectionValues,
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

    test("refreshes when the active source column changes data", () => {
        expect(
            shouldRefreshCategorySelectionValues(
                "data_changed",
                { columns: ["__tags"] },
                "__tags",
            ),
        ).toBe(true);

        expect(
            shouldRefreshCategorySelectionValues(
                "data_changed",
                { columns: ["other"] },
                "__tags",
            ),
        ).toBe(false);
    });

    test("refreshes when rows are added or the active source column is removed", () => {
        expect(
            shouldRefreshCategorySelectionValues(
                "data_added",
                10,
                "__tags",
            ),
        ).toBe(true);

        expect(
            shouldRefreshCategorySelectionValues(
                "column_removed",
                "__tags",
                "__tags",
            ),
        ).toBe(true);
    });
});
