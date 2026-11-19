import { describe, expect, test } from "vitest";
import {
    densityGridConfig,
    getDensityGridViewId,
    isDensityGridShown,
    shouldDrawInDensityGridCell,
    supportsDensityGridMode,
    withGridCrosshair,
} from "./densityGridUtils";

describe("densityGridUtils", () => {
    test("makes view ids safe and unique per cell", () => {
        expect(getDensityGridViewId("chart:#1", "marker/a", 2)).toBe("density-grid-chart__1-2-marker_a");
    });

    test("supports the deck scatter and SpatialData charts", () => {
        expect(supportsDensityGridMode("DeckContourScatter")).toBe(true);
        expect(supportsDensityGridMode("DeckDensity")).toBe(true);
        expect(supportsDensityGridMode("SpatialDataMdvRegionReact")).toBe(true);
        expect(supportsDensityGridMode("VivMdvRegionReact")).toBe(false);
        expect(supportsDensityGridMode("wgl_scatter_plot")).toBe(false);
        expect(supportsDensityGridMode(undefined)).toBe(false);
    });

    test("shows the grid only for a 2D chart with a grid by density fields", () => {
        const config = { type: "DeckContourScatter", dimension: "2d", grid: { by: "densityFields" } } as const;
        expect(isDensityGridShown(config)).toBe(true);
        expect(isDensityGridShown({ ...config, grid: null })).toBe(false);
        expect(isDensityGridShown({ ...config, grid: undefined })).toBe(false);
        expect(isDensityGridShown({ ...config, dimension: "3d" })).toBe(false);
        expect(isDensityGridShown({ ...config, dimension: undefined })).toBe(true);
        expect(isDensityGridShown({ ...config, type: "wgl_scatter_plot" })).toBe(false);
    });

    test("turns the grid on and off", () => {
        expect(densityGridConfig(true)).toEqual({ by: "densityFields" });
        expect(densityGridConfig(false)).toBeNull();
        // the grid's other options stay while it is on
        expect(densityGridConfig(true, { by: "densityFields", crosshair: true })).toEqual({
            by: "densityFields",
            crosshair: true,
        });
        expect(densityGridConfig(false, { by: "densityFields", crosshair: true })).toBeNull();
    });

    test("sets the crosshair with a new grid object, and only on a grid", () => {
        const grid = { by: "densityFields" } as const;
        const next = withGridCrosshair(grid, true);
        expect(next).toEqual({ by: "densityFields", crosshair: true });
        expect(next).not.toBe(grid);
        expect(withGridCrosshair(null, true)).toBeNull();
        expect(withGridCrosshair(undefined, true)).toBeUndefined();
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
