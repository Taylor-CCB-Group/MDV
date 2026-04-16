import { analyzeChartColumnImpact, describeColumnImpactReason } from "@/charts/columnRemovalUtils";
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
});
