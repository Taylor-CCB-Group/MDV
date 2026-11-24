import { describe, expect, it } from "vitest";
import { getMultiLineLayout, type MultiLineLayoutInput } from "@/react/components/multiline/multiLineLayout";

/** Every character is 0.6 of the font size wide. */
const measure = (text: string, fontSize: number) => text.length * fontSize * 0.6;

const base: MultiLineLayoutInput = {
    width: 300,
    height: 200,
    axis: { x: { size: 30, tickfont: 10 }, y: { size: 45, tickfont: 10 } },
    fitLabels: true,
    domain: [-1.5, 1.7],
    xLabel: "zscore_clipped_CD3",
    rowLabels: null,
    yLabel: "density",
};

const samples = ["Unchallenged-SAMPLE_6_ROI_1", "DNAonedose-Lesioned_SAMPLE_15_ROI_3", "Unchallenged-SAMPLE_6_ROI_2"];

describe("getMultiLineLayout", () => {
    it("fits the left margin to short row labels", () => {
        const layout = getMultiLineLayout({ ...base, rowLabels: ["0", "1", "10"] }, measure);
        expect(layout.rowLabels).toEqual(["0", "1", "10"]);
        // "10" is 12px wide, plus the tick and gap
        expect(layout.left).toBeLessThan(25);
    });

    it("shortens long row labels from the middle to leave most of the width for the plot", () => {
        const layout = getMultiLineLayout({ ...base, rowLabels: samples }, measure);
        expect(layout.plotWidth).toBeGreaterThanOrEqual(base.width * 0.6 - 15);
        for (const [i, label] of (layout.rowLabels ?? []).entries()) {
            expect(label).toContain("…");
            expect(label.startsWith(samples[i].slice(0, 3))).toBe(true);
            expect(label.endsWith(samples[i].slice(-3))).toBe(true);
            expect(measure(label, layout.rowFontSize)).toBeLessThanOrEqual(layout.left);
        }
    });

    it("keeps the legacy axis sizes when not fitting", () => {
        const layout = getMultiLineLayout({ ...base, fitLabels: false, rowLabels: samples }, measure);
        expect(layout.left).toBe(45);
        expect(layout.plotHeight).toBe(base.height - layout.top - 30);
        // what's left of 45px after the tick and gap
        for (const label of layout.rowLabels ?? []) expect(measure(label, 10)).toBeLessThanOrEqual(45);
    });

    it("shrinks row labels to the row height", () => {
        const many = Array.from({ length: 40 }, (_, i) => String(i));
        const layout = getMultiLineLayout({ ...base, rowLabels: many }, measure);
        expect(layout.rowFontSize).toBeLessThanOrEqual(layout.plotHeight / many.length);
    });

    it("uses fewer x ticks on a narrow plot, so their labels don't overlap", () => {
        const wide = getMultiLineLayout({ ...base, width: 800 }, measure);
        const narrow = getMultiLineLayout({ ...base, width: 120 }, measure);
        expect(narrow.xTicks.length).toBeLessThan(wide.xTicks.length);
        const spacing = narrow.plotWidth / (narrow.xTicks.length - 1);
        const widest = Math.max(...narrow.xTicks.map((t) => measure(narrow.xTickFormat(t), 10)));
        expect(widest).toBeLessThan(spacing);
    });

    it("thins to the outer ticks, then one, when even d3's fewest labels would overlap", () => {
        for (const width of [120, 60]) {
            const layout = getMultiLineLayout({ ...base, width, domain: [1000000, 1000001] }, measure);
            const { xTicks } = layout;
            expect(xTicks.length).toBeGreaterThan(0);
            const widest = Math.max(...xTicks.map((t) => measure(layout.xTickFormat(t), 10)));
            for (let i = 1; i < xTicks.length; i++) {
                const spacing = (layout.plotWidth * (xTicks[i] - xTicks[i - 1])) / 1;
                expect(widest).toBeLessThan(spacing);
            }
        }
    });

    it("shortens the x label to the plot width", () => {
        const layout = getMultiLineLayout({ ...base, width: 80 }, measure);
        expect(layout.xLabel).toContain("…");
        expect(measure(layout.xLabel, 10)).toBeLessThanOrEqual(layout.plotWidth);
    });
});
