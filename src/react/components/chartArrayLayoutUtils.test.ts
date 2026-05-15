import { describe, expect, test } from "vitest";
import {
    CHART_ARRAY_LAYOUT_DEFAULT_SIZING,
    chartArrayGridFitsViewport,
    estimateGridColumnCount,
    estimateSquareGridContentHeight,
    getChartArrayLayoutMode,
    getPreferredGridColumnCount,
} from "./chartArrayLayoutUtils";

const sizing = CHART_ARRAY_LAYOUT_DEFAULT_SIZING;

describe("getChartArrayLayoutMode", () => {
    test("uses single layout for zero or one cell", () => {
        expect(getChartArrayLayoutMode(0)).toBe("single");
        expect(getChartArrayLayoutMode(1)).toBe("single");
    });

    test("uses pair layout for two cells", () => {
        expect(getChartArrayLayoutMode(2)).toBe("pair");
    });

    test("uses grid layout for three or more cells", () => {
        expect(getChartArrayLayoutMode(3)).toBe("grid");
        expect(getChartArrayLayoutMode(12)).toBe("grid");
    });
});

describe("getPreferredGridColumnCount", () => {
    test("uses one row when every item can meet the minimum width", () => {
        expect(getPreferredGridColumnCount(4, 1200, sizing)).toBe(4);
    });

    test("wraps to fewer columns when a single-row share would be too narrow", () => {
        expect(getPreferredGridColumnCount(4, 800, sizing)).toBe(2);
        expect(estimateGridColumnCount(800, sizing)).toBe(3);
    });
});

describe("chart array grid viewport fit", () => {
    test("estimates max column count from min cell width", () => {
        expect(estimateGridColumnCount(500, sizing)).toBe(1);
        expect(estimateGridColumnCount(800, sizing)).toBe(3);
    });

    test("four cells in a wide short viewport need scrolling", () => {
        const innerWidth = 800;
        const contentHeight = estimateSquareGridContentHeight(4, innerWidth, sizing);
        expect(contentHeight).toBeGreaterThan(400);
        expect(chartArrayGridFitsViewport(4, 400, innerWidth, sizing)).toBe(false);
    });

    test("four cells in a tall viewport fit without scrolling", () => {
        const innerWidth = 800;
        const contentHeight = estimateSquareGridContentHeight(4, innerWidth, sizing);
        expect(chartArrayGridFitsViewport(4, contentHeight, innerWidth, sizing)).toBe(true);
        expect(chartArrayGridFitsViewport(4, contentHeight + 1, innerWidth, sizing)).toBe(true);
    });

    test("three cells in a typical panel height usually fit", () => {
        expect(chartArrayGridFitsViewport(3, 600, 800, sizing)).toBe(true);
    });
});
