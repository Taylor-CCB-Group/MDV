import { describe, expect, test, vi } from "vitest";
import type { BaseConfig } from "@/charts/BaseChart";
import type { MultiColumnQuery } from "@/links/link_utils";
import {
    analyzeChartColumnImpact,
    analyzeColumnRemoval,
    getReferencedFields,
} from "@/charts/columnRemovalUtils";

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
    test("classifies param usage as a blocking impact", () => {
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
                usagePaths: ["tooltip"],
            }),
        );
        expect(backgroundImpact).toEqual(
            expect.objectContaining({
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
                usage: "settings",
                usagePaths: ["tooltip_columns"],
            }),
        );
    });

    test("extracts referenced fields from nested field and fields objects", () => {
        const extracted = getReferencedFields({
            field: { field: "gene" },
            fields: [
                "age",
                { field: "cluster" },
                { fields: [{ field: "donor" }] },
                123,
            ],
            columnId: "sample_id",
        });

        expect(extracted).toEqual([
            "gene",
            "age",
            "cluster",
            "donor",
            "sample_id",
        ]);
    });

    test("detects nested config field references in analyzeChartColumnImpact", () => {
        const impact = analyzeChartColumnImpact(
            createChartConfig({
                title: "Nested",
                type: "custom",
                nested_config: { fields: [{ field: "age" }] },
            }),
            { configEntriesUsingColumns: ["nested_config"], name: "custom" },
            "age",
        );

        expect(impact).toEqual(
            expect.objectContaining({
                usage: "settings",
                usagePaths: ["nested_config"],
            }),
        );
    });

    test("collects impacts from current and saved views", async () => {
        const linkedQuery: MultiColumnQuery = {
            columns: [],
            fields: ["age"],
            initialize: async () => undefined,
        };
        const impact = await analyzeColumnRemoval({
            dataSourceName: "cells",
            columnName: "age",
            sourceChartId: "table-1",
            currentViewName: "Current View",
            allViewNames: ["Current View", "Other View"],
            currentCharts: [
                {
                    dataSourceName: "cells",
                    config: createChartConfig({
                        id: "table-1",
                        title: "",
                        type: "table_chart_react",
                        param: ["age"],
                    }),
                },
                {
                    dataSourceName: "cells",
                    config: createChartConfig({
                        id: "scatter-1",
                        title: "Scatter",
                        type: "wgl_scatter_plot_dev",
                        param: ["x", "y"],
                        color_by: "age",
                    }),
                },
                {
                    dataSourceName: "linked-ds",
                    config: createChartConfig({
                        id: "linked-1",
                        title: "Linked Scatter",
                        type: "wgl_scatter_plot_dev",
                        param: [linkedQuery, "y"],
                    }),
                },
            ],
            viewLoader: vi.fn(async () => ({
                initialCharts: {
                    cells: [
                        createChartConfig({
                            id: "saved-1",
                            title: "",
                            type: "selection_dialog",
                            param: ["age"],
                        }),
                    ],
                    linked: [
                        createChartConfig({
                            id: "saved-linked-1",
                            title: "Saved Linked",
                            type: "wgl_scatter_plot_dev",
                            param: [linkedQuery, "y"],
                        }),
                    ],
                },
            })),
        });

        expect(impact).toEqual({
            dataSourceName: "cells",
            columnName: "age",
            currentViewCharts: [
                expect.objectContaining({
                    chartId: "scatter-1",
                    usage: "settings",
                    usagePaths: ["color_by"],
                }),
            ],
            savedViews: [
                {
                    viewName: "Other View",
                    charts: [
                        expect.objectContaining({
                            chartId: "saved-1",
                            usage: "param",
                            usagePaths: ["param"],
                        }),
                    ],
                },
            ],
        });
    });

    test("fails closed when a saved view cannot be loaded", async () => {
        await expect(
            analyzeColumnRemoval({
                dataSourceName: "cells",
                columnName: "age",
                currentCharts: [],
                currentViewName: "Current View",
                allViewNames: ["Current View", "Broken View"],
                viewLoader: vi.fn(async () => {
                    throw new Error("load failed");
                }),
            }),
        ).rejects.toThrow(
            "Failed to check column usage in saved views: Broken View",
        );
    });
});
