import {
    analyzeChartColumnImpact,
    describeColumnImpactReason,
    getMissingColumnsForChartConfig,
    sanitizeChartConfigForRemovedColumns,
} from "@/charts/columnRemovalUtils";
import { describe, expect, test } from "vitest";

describe("analyzeChartColumnImpact", () => {
    test("deletes a chart when a single parameter is removed", () => {
        const impact = analyzeChartColumnImpact(
            {
                id: "chart-1",
                title: "Heat Map",
                type: "heat_map",
                size: [400, 300],
                legend: "",
                param: ["group", "expr1", "expr2"],
            } as any,
            {
                params: [
                    { type: "text", name: "group" },
                    { type: "_multi_column:number", name: "fields" },
                ],
            } as any,
            "group",
        );

        expect(impact?.action).toBe("delete");
        expect(impact?.reasons).toEqual([{ kind: "param.single" }]);
    });

    test("prunes a multi parameter and keeps the chart alive", () => {
        const impact = analyzeChartColumnImpact(
            {
                id: "chart-2",
                title: "Dot Plot",
                type: "dot_plot",
                size: [400, 300],
                legend: "",
                param: ["group", "expr1", "expr2"],
            } as any,
            {
                params: [
                    { type: "text", name: "group" },
                    { type: "_multi_column:number", name: "fields" },
                ],
            } as any,
            "expr1",
        );

        expect(impact?.action).toBe("update");
        expect(impact?.nextParam).toEqual(["group", "expr2"]);
        expect(impact?.paramSlotImpacts).toEqual([
            {
                paramIndex: 1,
                mode: "multi",
                action: "prune",
                remainingFields: 1,
            },
        ]);
    });

    test("deletes a chart when the last multi field is removed", () => {
        const impact = analyzeChartColumnImpact(
            {
                id: "chart-3",
                title: "Selection Dialog",
                type: "selection_dialog",
                size: [400, 300],
                legend: "",
                param: ["only_field"],
            } as any,
            {
                params: [{ type: "_multi_column:all", name: "fields" }],
            } as any,
            "only_field",
        );

        expect(impact?.action).toBe("delete");
        expect(impact?.paramSlotImpacts).toEqual([
            {
                paramIndex: 0,
                mode: "multi",
                action: "delete",
                remainingFields: 0,
            },
        ]);
    });

    test("resets color by, tooltip, and background filter without deleting", () => {
        const impact = analyzeChartColumnImpact(
            {
                id: "chart-4",
                title: "Scatter",
                type: "scatter",
                size: [400, 300],
                legend: "",
                param: ["x", "y"],
                color_by: "age",
                tooltip: {
                    show: true,
                    column: ["age", "name"],
                },
                background_filter: {
                    column: "age",
                    category: "A",
                },
            } as any,
            {
                params: [
                    { type: "number", name: "x" },
                    { type: "number", name: "y" },
                ],
            } as any,
            "age",
        );

        expect(impact?.action).toBe("update");
        expect(impact?.clearColorBy).toBe(true);
        expect(impact?.clearBackgroundFilter).toBe(true);
        expect(impact?.tooltipUpdate).toEqual({
            nextColumn: ["name"],
            disableTooltip: false,
        });
        expect(impact?.reasons.map(describeColumnImpactReason)).toEqual([
            "color by",
            "tooltip",
            "background filter",
        ]);
    });

    test("prunes array config entries and deletes scalar config entries", () => {
        const updateImpact = analyzeChartColumnImpact(
            {
                id: "chart-5",
                title: "Custom",
                type: "custom",
                size: [400, 300],
                legend: "",
                param: ["x"],
                densityFields: ["expr1", "expr2"],
            } as any,
            {
                params: [{ type: "number", name: "x" }],
                configEntriesUsingColumns: ["densityFields"],
            } as any,
            "expr1",
        );

        expect(updateImpact?.action).toBe("update");
        expect(updateImpact?.configEntryUpdates).toEqual({
            densityFields: ["expr2"],
        });

        const deleteImpact = analyzeChartColumnImpact(
            {
                id: "chart-6",
                title: "Custom",
                type: "custom",
                size: [400, 300],
                legend: "",
                param: ["x"],
                sortBy: "expr1",
            } as any,
            {
                params: [{ type: "number", name: "x" }],
                configEntriesUsingColumns: ["sortBy"],
            } as any,
            "expr1",
        );

        expect(deleteImpact?.action).toBe("delete");
        expect(deleteImpact?.reasons).toContainEqual({
            kind: "config_entry.single",
            entry: "sortBy",
        });
    });

    test("clears sort config entries without deleting table charts", () => {
        const impact = analyzeChartColumnImpact(
            {
                id: "chart-sort",
                title: "Table",
                type: "table_chart_react",
                size: [400, 300],
                legend: "",
                param: ["expr2"],
                sort: {
                    columnId: "expr1",
                    ascending: true,
                },
            } as any,
            {
                params: [{ type: "_multi_column:all", name: "fields" }],
                configEntriesUsingColumns: ["sort"],
            } as any,
            "expr1",
        );

        expect(impact?.action).toBe("update");
        expect(impact?.configEntryUpdates).toEqual({
            sort: null,
        });
        expect(impact?.reasons).toContainEqual({
            kind: "config_entry.single",
            entry: "sort",
        });
    });

    test("treats row summary image key removal as destructive", () => {
        const impact = analyzeChartColumnImpact(
            {
                id: "chart-7",
                title: "Row Summary",
                type: "row_summary_box",
                size: [400, 300],
                legend: "",
                param: ["field1", "image_key", "field2"],
                image: {
                    param: 1,
                },
            } as any,
            {
                params: [{ type: "_multi_column:all", name: "fields" }],
            } as any,
            "image_key",
        );

        expect(impact?.action).toBe("delete");
    });

    test("sanitizes selection dialog config when a removed field is pruned from a loaded view", () => {
        const sanitized = sanitizeChartConfigForRemovedColumns(
            {
                id: "chart-8",
                title: "Selection Dialog",
                type: "selection_dialog",
                size: [400, 300],
                legend: "",
                param: ["expr1", "expr2"],
                filters: {
                    expr1: null,
                    expr2: null,
                },
                order: {
                    expr1: 0,
                    expr2: 1,
                },
            } as any,
            {
                params: [{ type: "_multi_column:all", name: "fields" }],
            } as any,
            ["expr1"],
        );

        expect(sanitized).not.toBeNull();
        expect(sanitized?.param).toEqual(["expr2"]);
        expect((sanitized as any).filters).toEqual({ expr2: null });
        expect((sanitized as any).order).toEqual({ expr2: 0 });
    });

    test("sanitizes table config and removes stale sort/order/width state", () => {
        const sanitized = sanitizeChartConfigForRemovedColumns(
            {
                id: "chart-9",
                title: "Table",
                type: "table_chart_react",
                size: [400, 300],
                legend: "",
                param: ["expr1", "expr2"],
                order: {
                    expr1: 0,
                    expr2: 1,
                },
                sort: {
                    columnId: "expr1",
                    sortAsc: true,
                },
                column_widths: {
                    expr1: 120,
                    expr2: 180,
                },
            } as any,
            {
                params: [{ type: "_multi_column:all", name: "fields" }],
            } as any,
            ["expr1"],
        );

        expect(sanitized).not.toBeNull();
        expect(sanitized?.param).toEqual(["expr2"]);
        expect((sanitized as any).order).toEqual({ expr2: 0 });
        expect((sanitized as any).sort).toBeNull();
        expect((sanitized as any).column_widths).toEqual({ expr2: 180 });
    });

    test("shifts row summary image index when a preceding field is removed from a loaded view", () => {
        const sanitized = sanitizeChartConfigForRemovedColumns(
            {
                id: "chart-10",
                title: "Row Summary",
                type: "row_summary_box",
                size: [400, 300],
                legend: "",
                param: ["expr1", "image_key", "expr2"],
                image: {
                    param: 1,
                },
            } as any,
            {
                params: [{ type: "_multi_column:all", name: "fields" }],
            } as any,
            ["expr1"],
        );

        expect(sanitized).not.toBeNull();
        expect(sanitized?.param).toEqual(["image_key", "expr2"]);
        expect((sanitized as any).image.param).toBe(0);
    });

    test("finds missing columns referenced outside param", () => {
        const missingColumns = getMissingColumnsForChartConfig(
            {
                id: "chart-11",
                title: "Scatter",
                type: "scatter",
                size: [400, 300],
                legend: "",
                param: ["x", "y"],
                color_by: "removed_color",
                tooltip: {
                    show: true,
                    column: ["y", "removed_tooltip"],
                },
                background_filter: {
                    column: "removed_filter",
                },
            } as any,
            {
                params: [
                    { type: "number", name: "x" },
                    { type: "number", name: "y" },
                ],
            } as any,
            new Set(["x", "y"]),
        );

        expect(missingColumns.sort()).toEqual([
            "removed_color",
            "removed_filter",
            "removed_tooltip",
        ]);
    });

    test("finds missing sort columns referenced by tables", () => {
        const missingColumns = getMissingColumnsForChartConfig(
            {
                id: "chart-12",
                title: "Table",
                type: "table_chart_react",
                size: [400, 300],
                legend: "",
                param: ["expr2"],
                sort: {
                    columnId: "removed_sort",
                    ascending: true,
                },
            } as any,
            {
                params: [{ type: "_multi_column:all", name: "fields" }],
                configEntriesUsingColumns: ["sort"],
            } as any,
            new Set(["expr2"]),
        );

        expect(missingColumns).toEqual(["removed_sort"]);
    });
});
