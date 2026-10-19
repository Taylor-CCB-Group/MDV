import {
    buildFractionLegendSpec,
    getLegendBand,
    DOT_PLOT_LEGEND_MARGIN as m,
} from "@/react/components/dotplot/dotPlotLegends";
import { getFractionLegendSize } from "@/react/components/legend/FractionLegend";
import type { FractionLegendSpec } from "@/react/legend/fraction_legend/types";
import { describe, expect, it } from "vitest";

const color = { width: 195, height: 90 };
const fraction = { width: 260, height: 100 };

describe("getLegendBand", () => {
    it("puts the fraction legend beside the color legend when there's room", () => {
        const band = getLegendBand(1000, color, fraction);
        expect(band.color).toEqual([m, m]);
        expect(band.fraction).toEqual([m + color.width + m, m]);
        expect(band.height).toBe(m + fraction.height);
    });

    it("puts the fraction legend below the color legend when the chart is narrow", () => {
        const band = getLegendBand(300, color, fraction);
        expect(band.color).toEqual([m, m]);
        expect(band.fraction).toEqual([m, m + color.height + m]);
        expect(band.height).toBe(m + color.height + m + fraction.height);
    });

    it("only keeps space for the legends that are shown", () => {
        expect(getLegendBand(1000, color, null)).toEqual({ height: m + color.height, color: [m, m], fraction: null });
        expect(getLegendBand(1000, null, fraction)).toEqual({
            height: m + fraction.height,
            color: null,
            fraction: [m, m],
        });
        expect(getLegendBand(1000, null, null)).toEqual({ height: 0, color: null, fraction: null });
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
