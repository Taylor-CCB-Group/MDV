import { describe, expect, test } from "vitest";
import type { LoadedDataColumn } from "@/charts/charts";
import { buildCategoryHeatmapAggregation, pruneCategoryHeatmapAggregation } from "@/charts/CategoryHeatMap/categoryHeatmapUtils";

describe("categoryHeatmapUtils", () => {
    test("counts single-value categorical co-occurrences", () => {
        const xColumn = {
            datatype: "text",
            values: ["A", "B"],
            data: new Uint8Array([0, 1, 0]),
        } as LoadedDataColumn<"text">;
        const yColumn = {
            datatype: "text16",
            values: ["K", "L"],
            data: new Uint16Array([0, 1, 0]),
        } as LoadedDataColumn<"text16">;

        const result = buildCategoryHeatmapAggregation(
            xColumn,
            yColumn,
            new Uint32Array([0, 1, 2]),
        );

        expect(result.xLabels).toEqual(["A", "B"]);
        expect(result.yLabels).toEqual(["K", "L"]);
        expect(result.counts).toEqual([
            [2, 0],
            [0, 1],
        ]);
        expect(result.maxCount).toBe(2);
        expect(result.totalCells).toBe(4);
    });

    test("applies Cartesian pairing for multitext rows", () => {
        const xColumn = {
            datatype: "multitext",
            values: ["A", "B", "C"],
            data: new Uint8Array([0, 1, 2, 255]),
            stringLength: 2,
        } as unknown as LoadedDataColumn<"multitext">;
        const yColumn = {
            datatype: "multitext",
            values: ["K", "L", "M"],
            data: new Uint8Array([0, 1, 1, 2]),
            stringLength: 2,
        } as unknown as LoadedDataColumn<"multitext">;

        const result = buildCategoryHeatmapAggregation(
            xColumn,
            yColumn,
            new Uint32Array([0, 1]),
        );

        expect(result.xLabels).toEqual(["A", "B", "C"]);
        expect(result.yLabels).toEqual(["K", "L", "M"]);
        expect(result.counts).toEqual([
            [1, 1, 0],
            [1, 1, 1],
            [0, 0, 1],
        ]);
        expect(result.maxCount).toBe(1);
        expect(result.totalCells).toBe(9);
    });

    test("ignores rows where either side has no category", () => {
        const xColumn = {
            datatype: "text",
            values: ["A", ""],
            data: new Uint8Array([0, 1]),
        } as LoadedDataColumn<"text">;
        const yColumn = {
            datatype: "text",
            values: ["K", "L"],
            data: new Uint8Array([0, 1]),
        } as LoadedDataColumn<"text">;

        const result = buildCategoryHeatmapAggregation(
            xColumn,
            yColumn,
            new Uint32Array([0, 1]),
        );

        expect(result.xLabels).toEqual(["A"]);
        expect(result.yLabels).toEqual(["K"]);
        expect(result.counts).toEqual([[1]]);
        expect(result.maxCount).toBe(1);
    });

    test("prunes displayed categories without changing the original aggregation", () => {
        const xColumn = {
            datatype: "text",
            values: ["A", "B", "C"],
            data: new Uint8Array([0, 1, 2]),
        } as LoadedDataColumn<"text">;
        const yColumn = {
            datatype: "text",
            values: ["K", "L", "M"],
            data: new Uint8Array([0, 1, 2]),
        } as LoadedDataColumn<"text">;

        const aggregation = buildCategoryHeatmapAggregation(
            xColumn,
            yColumn,
            new Uint32Array([0, 1, 2]),
        );

        const pruned = pruneCategoryHeatmapAggregation(aggregation, {
            xCategories: ["A", "C"],
            yCategories: ["K"],
        });

        expect(aggregation.xLabels).toEqual(["A", "B", "C"]);
        expect(aggregation.yLabels).toEqual(["K", "L", "M"]);
        expect(pruned.xLabels).toEqual(["A", "C"]);
        expect(pruned.yLabels).toEqual(["K"]);
        expect(pruned.counts).toEqual([[1, 0]]);
        expect(pruned.maxCount).toBe(1);
    });
});
