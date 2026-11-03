import { getActiveChartScopeFilters } from "@/react/hooks";
import { describe, expect, test } from "vitest";

describe("getActiveChartScopeFilters", () => {
    test("includes a background filter with a column and categories", () => {
        const background_filter = { column: "cell_type", category: ["T"] };
        expect(getActiveChartScopeFilters({ background_filter })).toEqual([background_filter]);
    });

    test("ignores a background filter that is unset", () => {
        expect(getActiveChartScopeFilters({})).toEqual([]);
        expect(getActiveChartScopeFilters({ background_filter: null } as never)).toEqual([]);
        expect(getActiveChartScopeFilters({ background_filter: { column: null, category: ["T"] } } as never)).toEqual([]);
        expect(getActiveChartScopeFilters({ background_filter: { column: "__none__", category: ["T"] } })).toEqual([]);
        expect(getActiveChartScopeFilters({ background_filter: { column: "cell_type", category: [] } })).toEqual([]);
    });
});
