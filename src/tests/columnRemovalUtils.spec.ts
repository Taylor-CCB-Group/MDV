import { describe, expect, test } from "vitest";
import type { BaseConfig } from "@/charts/BaseChart";
import { analyzeChartColumnImpact } from "@/charts/columnRemovalUtils";

type TestChartConfig = BaseConfig & Record<string, unknown>;

function createChartConfig(overrides: Partial<TestChartConfig> = {}): TestChartConfig {
    return {
        id: "chart-1",
        title: "Chart",
        type: "scatter_plot",
        param: ["x", "y"] as BaseConfig["param"],
        size: [400, 300],
        legend: "none",
        ...overrides,
    };
}

describe("columnRemovalUtils", () => {
    test("classifies param usage as chart removal", () => {
        const impact = analyzeChartColumnImpact(
            createChartConfig({
                title: "Scatter",
                param: ["x", "age"],
            }),
            undefined,
            "age",
        );

        expect(impact).toEqual(
            expect.objectContaining({
                action: "remove_chart",
                usage: "param",
                chartTitle: "Scatter",
                usagePaths: ["param"],
            }),
        );
    });

    test("classifies color_by usage as blocking delete", () => {
        const impact = analyzeChartColumnImpact(
            createChartConfig({
                title: "Scatter",
                color_by: "age",
            }),
            undefined,
            "age",
        );

        expect(impact).toEqual(
            expect.objectContaining({
                action: "block_delete",
                usage: "settings",
                usagePaths: ["color_by"],
            }),
        );
    });

    test("classifies tooltip and background filter usage as blocking delete", () => {
        const tooltipImpact = analyzeChartColumnImpact(
            createChartConfig({
                title: "Scatter",
                tooltip: { column: "age" },
            }),
            undefined,
            "age",
        );
        const backgroundImpact = analyzeChartColumnImpact(
            createChartConfig({
                id: "chart-2",
                title: "Heatmap",
                type: "heat_map",
                background_filter: { column: "age" },
            }),
            undefined,
            "age",
        );

        expect(tooltipImpact).toEqual(
            expect.objectContaining({
                action: "block_delete",
                usagePaths: ["tooltip"],
            }),
        );
        expect(backgroundImpact).toEqual(
            expect.objectContaining({
                action: "block_delete",
                usagePaths: ["background_filter"],
            }),
        );
    });

    test("classifies configEntriesUsingColumns usage as blocking delete", () => {
        const impact = analyzeChartColumnImpact(
            createChartConfig({
                title: "Custom",
                type: "custom",
                tooltip_columns: [{ field: "age" }],
            }),
            { configEntriesUsingColumns: ["tooltip_columns"], name: "custom" },
            "age",
        );

        expect(impact).toEqual(
            expect.objectContaining({
                action: "block_delete",
                usage: "settings",
                usagePaths: ["tooltip_columns"],
            }),
        );
    });

});
