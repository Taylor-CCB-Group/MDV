import {
    DOT_PLOT_MAX_RADIUS,
    DOT_PLOT_RIGHT_MARGIN,
    DOT_PLOT_TOP_MARGIN,
    getDotPlotLayout,
    getDotRadius,
    getVisibleRows,
    getXLabelLayout,
} from "@/react/components/dotplot/dotPlotLayout";
import { MISSING_VALUE_COLOR, binnedColorScale } from "@/react/components/dotplot/dotPlotPalettes";
import type { DotPlotGroup } from "@/react/components/dotplot/dotPlotStats";
import { describe, expect, it } from "vitest";

const group = (count: number): DotPlotGroup => ({ count, cells: [] });

describe("getXLabelLayout", () => {
    // every character 6px wide at size 10, scaling with font size
    const measure = (text: string, fontSize: number) => text.length * fontSize * 0.6;
    const axis = {
        size: 60,
        tickfont: 10,
        rotate_labels: true,
        fit_labels: true,
        truncate_labels: true,
        thin_labels: false,
    };
    const space = { plotLeft: 100, plotWidth: 300, chartHeight: 400 };

    it("shortens long labels to the set axis height when auto-fit is off, if asked", () => {
        const long = "ENSG00000228253_readthrough_long_noncoding";
        const off = { ...axis, fit_labels: false };
        const shortened = getXLabelLayout(off, [long, long], space, measure);
        expect(shortened.margin).toBe(60);
        expect(shortened.labels[0]).toContain("…");
        expect(getXLabelLayout({ ...off, truncate_labels: false }, [long, long], space, measure).labels).toEqual([
            long,
            long,
        ]);
    });

    it("uses the axis settings as they are when auto-fit is off", () => {
        const labels = ["GZMB", "NKG7"];
        expect(getXLabelLayout({ ...axis, fit_labels: false }, labels, space, measure)).toEqual({
            angle: 45,
            fontSize: 10,
            margin: 60,
            labels,
        });
        expect(
            getXLabelLayout({ ...axis, fit_labels: false, rotate_labels: false }, labels, space, measure).angle,
        ).toBe(0);
    });

    it("lays short labels flat, with only as much height as they need", () => {
        const placement = getXLabelLayout(axis, ["GZMB", "NKG7", "CD3E"], space, measure);
        expect(placement.angle).toBe(0);
        expect(placement.margin).toBeLessThan(axis.size);
        expect(placement.labels).toEqual(["GZMB", "NKG7", "CD3E"]);
    });

    it("rotates and shortens long labels when they don't fit", () => {
        const long = "ENSG00000228253_readthrough_long_noncoding";
        const placement = getXLabelLayout(axis, [long, long, long, long, long, long], space, measure);
        expect(placement.angle).toBe(45);
        expect(placement.labels[0]).toContain("…");
        expect(placement.labels[0].length).toBeLessThan(long.length);
    });

    it("keeps full labels with shortening off", () => {
        const long = "ENSG00000228253_readthrough_long_noncoding";
        const placement = getXLabelLayout({ ...axis, truncate_labels: false }, [long, long, long], space, measure);
        expect(placement.labels).toEqual([long, long, long]);
    });

    it("leaves some fields unlabelled when hiding overlapping labels in a crowded axis", () => {
        const labels = Array.from({ length: 60 }, (_, i) => `GENE${i}`);
        const thinned = getXLabelLayout({ ...axis, thin_labels: true }, labels, space, measure);
        expect(thinned.labels.some((l) => l === "")).toBe(true);
        const all = getXLabelLayout(axis, labels, space, measure);
        expect(all.labels.every((l) => l !== "")).toBe(true);
    });

    it("falls back to the axis settings when there's no room to fit", () => {
        expect(getXLabelLayout(axis, ["GZMB"], { ...space, plotWidth: 0 }, measure).margin).toBe(60);
    });
});

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

    it("drops the panel gaps when there's no room, so panels stay inside the chart", () => {
        // 10px left for three panels and two 20px gaps
        const layout = getDotPlotLayout({ ...input, width: 120, panelCount: 3 });
        expect(layout.panelWidth).toBeCloseTo(10 / 3);
        for (const left of layout.panelLefts) {
            expect(left + layout.panelWidth).toBeLessThanOrEqual(120 - DOT_PLOT_RIGHT_MARGIN + 1e-9);
        }
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
