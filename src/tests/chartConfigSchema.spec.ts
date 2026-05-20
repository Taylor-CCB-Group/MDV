import { describe, expect, test } from "vitest";
import { getChartConfigSchema } from "@/charts/schemas/ChartConfigRegistry";
import { safeValidateChartConfig } from "@/charts/schemas/ChartConfigSchema";

describe("chart config schema registration", () => {
    test("validates text box configs for runtime and legacy aliases", () => {
        const runtimeConfig = {
            id: "text-box-runtime",
            title: "Text Box Runtime",
            type: "text_box_chart",
            text: "## Section\nRuntime config",
            param: [],
            size: [420, 240],
        };

        const legacyAliasConfig = {
            ...runtimeConfig,
            id: "text-box-legacy",
            title: "Text Box Legacy",
            type: "text_box",
        };

        expect(getChartConfigSchema("text_box_chart")).toBeDefined();
        expect(getChartConfigSchema("text_box")).toBeDefined();
        expect(safeValidateChartConfig(runtimeConfig)).toMatchObject({
            type: "text_box_chart",
            text: "## Section\nRuntime config",
        });
        expect(safeValidateChartConfig(legacyAliasConfig)).toMatchObject({
            type: "text_box",
            text: "## Section\nRuntime config",
        });
    });

    test("validates the new deck contour scatter chart config", () => {
        const config = {
            id: "deck-contour-1",
            title: "Contours",
            type: "DeckContourScatter",
            param: ["x", "y"],
            size: [640, 480],
            axis: {
                x: { size: 20, tickfont: 10, rotate_labels: false },
                y: { size: 40, tickfont: 10, rotate_labels: false },
            },
            viewState: {
                target: [0, 0, 0],
                zoom: 0,
                minZoom: -50,
            },
            tooltip: {
                show: true,
                column: ["label"],
            },
            contour_fill: false,
            contour_bandwidth: 2,
            contour_intensity: 0.8,
            contour_opacity: 0.5,
            contour_fillThreshold: 2,
            contourParameter: "cluster",
            category1: ["A"],
            densityFields: ["marker_a", "marker_b"],
            field_legend: {
                display: true,
            },
            category_filters: [
                {
                    column: "sample",
                    category: ["all"],
                },
            ],
            point_shape: "gaussian",
            on_filter: "grey",
            zoom_on_filter: true,
        };

        expect(getChartConfigSchema("DeckContourScatter")).toBeDefined();
        expect(safeValidateChartConfig(config)).toMatchObject({
            type: "DeckContourScatter",
            contourParameter: "cluster",
            densityFields: ["marker_a", "marker_b"],
        });
    });

    test("validates the legacy DeckDensity alias against the same config family", () => {
        const config = {
            id: "deck-density-1",
            title: "Legacy density",
            type: "DeckDensity",
            param: ["x", "y"],
            size: [320, 240],
            viewState: {
                target: [0, 0, 0],
                zoom: 1,
            },
            tooltip: {
                show: false,
            },
            contour_fill: true,
            contour_bandwidth: 1,
            contour_intensity: 0.5,
            contour_opacity: 0.3,
            contour_fillThreshold: 1.5,
            field_legend: {
                display: false,
            },
        };

        expect(getChartConfigSchema("DeckDensity")).toBeDefined();
        expect(safeValidateChartConfig(config)).toMatchObject({
            type: "DeckDensity",
            contour_fill: true,
        });
    });

    test("validates DeckSplatter against the shared contour and density field config sections", () => {
        const config = {
            id: "deck-splatter-1",
            title: "Splatter",
            type: "DeckSplatter",
            param: ["x", "y"],
            size: [500, 400],
            category: "cell_type",
            densityFields: ["marker_a", "marker_b"],
            contour_fill: true,
            contour_bandwidth: 1.5,
            contour_intensity: 0.75,
            contour_opacity: 0.4,
            contour_fillThreshold: 2,
        };

        expect(getChartConfigSchema("DeckSplatter")).toBeDefined();
        expect(safeValidateChartConfig(config)).toMatchObject({
            type: "DeckSplatter",
            category: "cell_type",
            densityFields: ["marker_a", "marker_b"],
            contour_bandwidth: 1.5,
        });
    });

    test("validates category heatmap config", () => {
        const config = {
            id: "category-heatmap-1",
            title: "Category Heatmap",
            type: "category_heatmap",
            param: ["cluster", "sample_type"],
            size: [640, 400],
        };

        expect(getChartConfigSchema("category_heatmap")).toBeDefined();
        expect(safeValidateChartConfig(config)).toMatchObject({
            type: "category_heatmap",
            param: ["cluster", "sample_type"],
        });
    });
});
