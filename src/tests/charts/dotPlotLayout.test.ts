import {
    DOT_PLOT_MAX_RADIUS,
    DOT_PLOT_RIGHT_MARGIN,
    DOT_PLOT_TOP_MARGIN,
    getDotPlotLayout,
    getDotRadius,
    getVisibleRows,
} from "@/react/components/dotplot/dotPlotLayout";
import { MISSING_VALUE_COLOR, binnedColorScale } from "@/react/components/dotplot/dotPlotPalettes";
import type { DotPlotGroup } from "@/react/components/dotplot/dotPlotStats";
import { describe, expect, it } from "vitest";

const group = (count: number): DotPlotGroup => ({ count, cells: [] });

describe("getVisibleRows", () => {
    const values = ["T-cell", "b-cell", "NK-cell", "B10", "B2"];

    it("drops rows with no data in any panel, keeping data order", () => {
        const panels = [[group(3), group(0), group(1), group(0), group(2)]];
        expect(getVisibleRows(panels, values, "data")).toEqual([0, 2, 4]);
    });

    it("keeps a row that has data in only one of the panels", () => {
        const panels = [
            [group(0), group(0), group(1), group(0), group(0)],
            [group(2), group(0), group(0), group(0), group(0)],
        ];
        expect(getVisibleRows(panels, values, "data")).toEqual([0, 2]);
    });

    it("sorts alphabetically, ignoring case and with numbers in numeric order", () => {
        const panels = [values.map(() => group(1))];
        const order = getVisibleRows(panels, values, "alphabetical").map((r) => values[r]);
        expect(order).toEqual(["b-cell", "B2", "B10", "NK-cell", "T-cell"]);
    });
});

describe("getDotPlotLayout", () => {
    const input = {
        width: 600,
        height: 400,
        yAxisWidth: 100,
        xAxisHeight: 60,
        titleHeight: 0,
        panelGap: 20,
        panelCount: 1,
        fieldCount: 5,
        rowCount: 4,
    };

    it("fills the space left of the axes with one panel", () => {
        const layout = getDotPlotLayout(input);
        expect(layout.top).toBe(DOT_PLOT_TOP_MARGIN);
        expect(layout.panelLefts).toEqual([100]);
        expect(layout.panelWidth).toBe(600 - 100 - DOT_PLOT_RIGHT_MARGIN);
        expect(layout.plotHeight).toBe(400 - DOT_PLOT_TOP_MARGIN - 60);
        expect(layout.cellWidth).toBe(layout.panelWidth / 5);
        expect(layout.cellHeight).toBe(layout.plotHeight / 4);
    });

    it("splits the width evenly between panels, with gaps between them", () => {
        const layout = getDotPlotLayout({ ...input, panelCount: 3, titleHeight: 20 });
        const width = (600 - 100 - DOT_PLOT_RIGHT_MARGIN - 2 * 20) / 3;
        expect(layout.panelWidth).toBeCloseTo(width);
        expect(layout.panelLefts[0]).toBe(100);
        expect(layout.panelLefts[1]).toBeCloseTo(100 + width + 20);
        expect(layout.panelLefts[2]).toBeCloseTo(100 + 2 * (width + 20));
        expect(layout.top).toBe(DOT_PLOT_TOP_MARGIN + 20);
    });

    it("fits dots in the smaller cell side and caps the radius", () => {
        expect(getDotPlotLayout({ ...input, fieldCount: 50 }).maxRadius).toBeCloseTo(490 / 50 / 2);
        expect(getDotPlotLayout({ ...input, fieldCount: 1, rowCount: 1 }).maxRadius).toBe(DOT_PLOT_MAX_RADIUS);
    });

    it("never goes negative when the chart is too small", () => {
        const layout = getDotPlotLayout({ ...input, width: 50, height: 30 });
        expect(layout.panelWidth).toBe(0);
        expect(layout.plotHeight).toBe(0);
        expect(layout.maxRadius).toBe(0);
    });
});

describe("getDotRadius", () => {
    it("scales dot area with the fraction", () => {
        expect(getDotRadius(100, 100, 20)).toBe(20);
        expect(getDotRadius(25, 100, 20)).toBe(10);
        expect(getDotRadius(0, 100, 20)).toBe(0);
        expect(getDotRadius(50, 0, 20)).toBe(0);
    });
});

describe("binnedColorScale", () => {
    // as from DataStore.getColumnColors with 4 bins: one colour per bin edge
    const colors = ["c0", "c1", "c2", "c3", "c4"];
    const scale = binnedColorScale(colors, 0, 4);

    it("bins values evenly between min and max, clamping outside", () => {
        expect(scale(0)).toBe("c0");
        expect(scale(1.5)).toBe("c1");
        expect(scale(3.99)).toBe("c3");
        expect(scale(4)).toBe("c4");
        expect(scale(-1)).toBe("c0");
        expect(scale(10)).toBe("c4");
    });

    it("uses the missing colour for NaN and the first colour when min equals max", () => {
        expect(scale(Number.NaN)).toBe(MISSING_VALUE_COLOR);
        expect(binnedColorScale(colors, 2, 2)(2)).toBe("c0");
    });

    it("accepts RGB arrays", () => {
        expect(binnedColorScale([[1, 2, 3]], 0, 1)(0)).toBe("rgb(1,2,3)");
    });
});
