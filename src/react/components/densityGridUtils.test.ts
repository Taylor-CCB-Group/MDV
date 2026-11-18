import { describe, expect, test } from "vitest";
import {
    getDensityGridViewId,
    isDensityGridShown,
    shouldDrawInDensityGridCell,
    supportsDensityGridMode,
} from "./densityGridUtils";

describe("densityGridUtils", () => {
    test("makes view ids safe and unique per cell", () => {
        expect(getDensityGridViewId("chart:#1", "marker/a", 2)).toBe("density-grid-chart__1-2-marker_a");
    });

    test("supports the deck scatter types with density fields", () => {
        expect(supportsDensityGridMode("DeckContourScatter")).toBe(true);
        expect(supportsDensityGridMode("DeckDensity")).toBe(true);
        expect(supportsDensityGridMode("VivMdvRegionReact")).toBe(false);
        expect(supportsDensityGridMode("wgl_scatter_plot")).toBe(false);
        expect(supportsDensityGridMode(undefined)).toBe(false);
    });

    test("shows the grid only for a 2D chart set to grid", () => {
        const config = { type: "DeckContourScatter", dimension: "2d", density_mode: "grid" } as const;
        expect(isDensityGridShown(config)).toBe(true);
        expect(isDensityGridShown({ ...config, density_mode: "overlay" })).toBe(false);
        expect(isDensityGridShown({ ...config, density_mode: undefined })).toBe(false);
        expect(isDensityGridShown({ ...config, dimension: "3d" })).toBe(false);
        expect(isDensityGridShown({ ...config, type: "wgl_scatter_plot" })).toBe(false);
    });

    test("draws a cell's layers and their sublayers only in that cell", () => {
        const cellA = getDensityGridViewId("chart-1", "field_a", 0);
        const cellB = getDensityGridViewId("chart-1", "field_b", 1);
        expect(shouldDrawInDensityGridCell(`${cellA}-density`, cellA)).toBe(true);
        expect(shouldDrawInDensityGridCell(`${cellA}-density-heatmap`, cellA)).toBe(true);
        expect(shouldDrawInDensityGridCell(`${cellA}-density`, cellB)).toBe(false);
    });

    test("does not confuse cells whose indices share a prefix", () => {
        const cell1 = getDensityGridViewId("c", "x", 1);
        const cell10 = getDensityGridViewId("c", "x", 10);
        expect(shouldDrawInDensityGridCell(`${cell10}-density`, cell1)).toBe(false);
        expect(shouldDrawInDensityGridCell(`${cell1}-density`, cell10)).toBe(false);
    });

    test("draws every other layer in every cell", () => {
        const cell = getDensityGridViewId("chart-1", "field_a", 0);
        expect(shouldDrawInDensityGridCell("scatter_-#chart-1detail-react#", cell)).toBe(true);
        expect(shouldDrawInDensityGridCell("selection_-#chart-1detail-react#-guides", cell)).toBe(true);
    });
});
