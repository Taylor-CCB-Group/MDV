import { describe, expect, test } from "vitest";
import { analyzeChartColumnImpact } from "@/charts/columnRemovalUtils";

function createChartConfig(overrides: Record<string, unknown> = {}) {
    return {
        id: "chart-1",
        title: "Chart",
        type: "scatter_plot",
        param: ["x", "y"],
        size: [400, 300] as [number, number],
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
            }) as any,
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
            }) as any,
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
            }) as any,
            undefined,
            "age",
        );
        const backgroundImpact = analyzeChartColumnImpact(
            createChartConfig({
                id: "chart-2",
                title: "Heatmap",
                type: "heat_map",
                background_filter: { column: "age" },
            }) as any,
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
            }) as any,
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
