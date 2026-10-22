import {
    MAX_LEGEND_BAND_FRACTION,
    buildFractionLegendSpec,
    getColorLegendBox,
    getFractionLegendBox,
    getLegendLayout,
    DOT_PLOT_LEGEND_MARGIN as m,
    pickEvenly,
} from "@/react/components/dotplot/dotPlotLegends";
import { getFractionLegendSize } from "@/react/components/legend/FractionLegend";
import type { FractionLegendSpec } from "@/react/legend/fraction_legend/types";
import { describe, expect, it } from "vitest";

const colorRange: [number, number] = [0, 90];
const input = { width: 1200, height: 800, colorRange, fracMax: 100, maxRadius: 25 };
/** Width of a one-row band at full size: colour legend, then the fraction legend beside it. */
const fullRowWidth =
    m + getColorLegendBox(colorRange).width + m + getFractionLegendBox(buildFractionLegendSpec(100, 25)).width + m;

describe("getLegendLayout", () => {
    it("puts both legends at full size in one row when there's room", () => {
        const layout = getLegendLayout(input);
        expect(layout.fractionItems).toBe(Number.POSITIVE_INFINITY);
        expect(layout.color).toEqual([m, m]);
        expect(layout.fraction?.[1]).toBe(m);
        expect(layout.fraction?.[0]).toBeGreaterThan(m);
        expect(layout.height).toBeGreaterThan(0);
    });

    it("shows fewer fraction dots rather than stacking when the row is a little too wide", () => {
        const layout = getLegendLayout({ ...input, width: fullRowWidth - 1 });
        expect(layout.fractionItems).toBeLessThan(Number.POSITIVE_INFINITY);
        expect(layout.fraction?.[1]).toBe(m);
        expect(layout.height).toBeGreaterThan(0);
    });

    it("never stacks the legends in the band", () => {
        for (const width of [300, 400, 500, 600, 800]) {
            const layout = getLegendLayout({ ...input, width });
            if (layout.height > 0) {
                expect(layout.fraction?.[1]).toBe(m);
            }
        }
    });

    it("drops the band and starts the legends over the plot when even a compact row won't fit", () => {
        const layout = getLegendLayout({ ...input, width: 250 });
        expect(layout.height).toBe(0);
        expect(layout.color).toEqual([m, m]);
        // no room beside the colour legend, so the fraction legend starts below it
        expect(layout.fraction?.[0]).toBe(m);
        expect(layout.fraction?.[1]).toBeGreaterThan(m);
    });

    it("never lets the band take more than a third of the height", () => {
        const layout = getLegendLayout({ ...input, height: 150 });
        expect(layout.height).toBe(0);
        for (const height of [300, 400, 800]) {
            expect(getLegendLayout({ ...input, height }).height).toBeLessThanOrEqual(height * MAX_LEGEND_BAND_FRACTION);
        }
    });

    it("only lays out the legends that are shown", () => {
        const colorOnly = getLegendLayout({ ...input, fracMax: null });
        expect(colorOnly.fraction).toBeNull();
        expect(colorOnly.color).toEqual([m, m]);
        const fractionOnly = getLegendLayout({ ...input, colorRange: null });
        expect(fractionOnly.color).toBeNull();
        expect(fractionOnly.fraction).toEqual([m, m]);
        const none = getLegendLayout({ ...input, colorRange: null, fracMax: null });
        expect(none).toMatchObject({ height: 0, color: null, fraction: null });
    });
});

describe("pickEvenly", () => {
    it("keeps the first and last, spread evenly", () => {
        expect(pickEvenly([20, 40, 60, 80, 100], 3)).toEqual([20, 60, 100]);
        expect(pickEvenly([20, 40, 60, 80, 100], 2)).toEqual([20, 100]);
    });

    it("keeps everything when asking for as many or more", () => {
        expect(pickEvenly([1, 2], 5)).toEqual([1, 2]);
        expect(pickEvenly([1, 2, 3], Number.POSITIVE_INFINITY)).toEqual([1, 2, 3]);
    });

    it("keeps the largest when asking for one", () => {
        expect(pickEvenly([20, 40, 60], 1)).toEqual([60]);
    });
});

describe("buildFractionLegendSpec", () => {
    it("makes a horizontal legend with dots the same size as the panels' dots", () => {
        const spec = buildFractionLegendSpec(100, 20);
        expect(spec.orientation).toBe("horizontal");
        expect(spec.items.map((i) => i.label)).toEqual(["20", "40", "60", "80", "100"]);
        expect(spec.items[spec.items.length - 1].radius).toBe(20);
        expect(spec.items[0].radius).toBeCloseTo(20 * Math.sqrt(0.2));
    });

    it("shows fewer dots when asked, keeping the smallest and largest", () => {
        expect(buildFractionLegendSpec(100, 20, 3).items.map((i) => i.label)).toEqual(["20", "60", "100"]);
    });

    it("has no items when nothing is expressed", () => {
        expect(buildFractionLegendSpec(0, 20).items).toEqual([]);
    });
});

describe("getFractionLegendSize", () => {
    const items = [
        { key: "20", label: "20", radius: 5 },
        { key: "80", label: "80", radius: 10 },
    ];

    it("keeps the vertical legend's fixed width", () => {
        const spec: FractionLegendSpec = { label: "fraction", items };
        expect(getFractionLegendSize(spec).width).toBe(100);
    });

    it("sizes a horizontal legend to one row of dots", () => {
        const vertical = getFractionLegendSize({ label: "% expressing", items });
        const horizontal = getFractionLegendSize({ label: "% expressing", items, orientation: "horizontal" });
        expect(horizontal.height).toBeLessThan(vertical.height);
        // two dots and their labels side by side
        expect(horizontal.width).toBeGreaterThan(2 * 5 + 2 * 10);
    });

    it("respects an explicit width and maximum height", () => {
        const size = getFractionLegendSize({ label: "x", items, orientation: "horizontal", width: 50, maxHeight: 30 });
        expect(size).toEqual({ width: 50, height: 30 });
    });
});
