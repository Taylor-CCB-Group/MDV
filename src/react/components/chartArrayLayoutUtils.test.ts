import { describe, expect, test } from "vitest";
import {
    CHART_ARRAY_LAYOUT_DEFAULT_SIZING,
    chartArrayCellSize,
    chooseChartArrayGrid,
    estimateGridColumnCount,
    getChartArrayLayoutMode,
    getPreferredGridColumnCount,
} from "./chartArrayLayoutUtils";

const sizing = CHART_ARRAY_LAYOUT_DEFAULT_SIZING;

describe("getChartArrayLayoutMode", () => {
    test("uses single layout for zero or one cell", () => {
        expect(getChartArrayLayoutMode(0)).toBe("single");
        expect(getChartArrayLayoutMode(1)).toBe("single");
    });

    test("uses grid layout for two or more cells", () => {
        expect(getChartArrayLayoutMode(2)).toBe("grid");
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

    test("estimates max column count from min cell width", () => {
        expect(estimateGridColumnCount(500, sizing)).toBe(1);
        expect(estimateGridColumnCount(800, sizing)).toBe(3);
    });
});

describe("chooseChartArrayGrid", () => {
    test("a wide, short panel puts the cells in one row", () => {
        expect(chooseChartArrayGrid(4, 1400, 320, sizing)).toEqual({ columns: 4, fillsViewport: true });
    });

    test("a tall panel stacks two cells rather than splitting it into two tall columns", () => {
        expect(chooseChartArrayGrid(2, 600, 1000, sizing)).toEqual({ columns: 1, fillsViewport: true });
        expect(chooseChartArrayGrid(2, 1000, 600, sizing)).toEqual({ columns: 2, fillsViewport: true });
    });

    test("a square panel makes a grid, not a row of tall columns", () => {
        // one row would be 287px wide and 876px tall
        expect(chooseChartArrayGrid(3, 876, 876, sizing)).toEqual({ columns: 2, fillsViewport: true });
        expect(chooseChartArrayGrid(4, 1200, 1200, sizing)).toEqual({ columns: 2, fillsViewport: true });
    });

    test("prefers no empty slots when the cells come out about the same size", () => {
        // 2 columns: 596 x 392; 3 columns: 394 x 392 with two empty slots
        expect(chooseChartArrayGrid(4, 1200, 792, sizing).columns).toBe(2);
    });

    test("scrolls square cells when they cannot all fit at the minimum size", () => {
        expect(chooseChartArrayGrid(40, 1200, 800, sizing)).toEqual({
            columns: getPreferredGridColumnCount(40, 1200, sizing),
            fillsViewport: false,
        });
        expect(chooseChartArrayGrid(4, 800, 400, sizing).fillsViewport).toBe(false);
    });

    test("no cells, or no width yet, is one column that fills", () => {
        expect(chooseChartArrayGrid(0, 800, 800, sizing)).toEqual({ columns: 1, fillsViewport: true });
        expect(chooseChartArrayGrid(4, 0, 800, sizing)).toEqual({ columns: 1, fillsViewport: true });
    });
});

describe("chartArrayCellSize", () => {
    test("rounds a filling cell down to whole device pixels each way", () => {
        // (1117 - 8) / 2 = 554.5 CSS px, 1109 device px at 2x; (563 - 8) / 2 = 277.5
        const size = chartArrayCellSize(4, { columns: 2, fillsViewport: true }, 1117, 563, 2, sizing);
        expect(size).toEqual({ width: 554.5, height: 277.5 });
        const odd = chartArrayCellSize(3, { columns: 3, fillsViewport: true }, 1000, 600, 1.5, sizing);
        expect(odd.width * 1.5).toBe(Math.floor(((1000 - 16) / 3) * 1.5));
        expect(Number.isInteger(odd.height * 1.5)).toBe(true);
    });

    test("makes a scrolling cell square", () => {
        const size = chartArrayCellSize(9, { columns: 2, fillsViewport: false }, 601, 300, 1, sizing);
        expect(size).toEqual({ width: 296, height: 296 });
    });

    test("treats a missing pixel ratio as 1", () => {
        expect(chartArrayCellSize(2, { columns: 2, fillsViewport: true }, 601, 400, 0, sizing)).toEqual({
            width: 296,
            height: 400,
        });
    });
});
