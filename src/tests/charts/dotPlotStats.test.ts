import { type DotPlotStatsRequest, computeDotPlotStats } from "@/react/components/dotplot/dotPlotStats";
import { describe, expect, it } from "vitest";

// six data rows: three of row category 0, three of row category 1
const base: DotPlotStatsRequest = {
    rowCodes: new Uint8Array([0, 0, 0, 1, 1, 1]),
    rowCount: 2,
    fields: [new Float32Array([2, 0, 4, 1, Number.NaN, 0]), new Int32Array([0, 0, 0, 3, 3, 3])],
};

const cell = (stats: ReturnType<typeof computeDotPlotStats>, panel: number, row: number, field: number) =>
    stats.panels[panel][row].cells[field];

describe("computeDotPlotStats", () => {
    it("averages over every row in the group, counting unexpressed rows as zero", () => {
        const stats = computeDotPlotStats(base);
        expect(stats.panels).toHaveLength(1);
        expect(stats.panels[0].map((g) => g.count)).toEqual([3, 3]);
        expect(cell(stats, 0, 0, 0).mean).toBeCloseTo(2);
        expect(cell(stats, 0, 0, 0).frac).toBeCloseTo(200 / 3);
        expect(cell(stats, 0, 0, 1)).toEqual({ mean: 0, frac: 0 });
        expect(cell(stats, 0, 1, 1)).toEqual({ mean: 3, frac: 100 });
        expect(stats.meanRange).toEqual([0, 3]);
        expect(stats.fracMax).toBe(100);
    });

    it("treats NaN as unexpressed but keeps the row in the group count", () => {
        const stats = computeDotPlotStats(base);
        expect(cell(stats, 0, 1, 0).mean).toBeCloseTo(1 / 3);
        expect(cell(stats, 0, 1, 0).frac).toBeCloseTo(100 / 3);
    });

    it("only counts values above the threshold", () => {
        const stats = computeDotPlotStats({ ...base, threshold: 1.5 });
        expect(cell(stats, 0, 0, 0).mean).toBeCloseTo(2);
        expect(cell(stats, 0, 1, 0)).toEqual({ mean: 0, frac: 0 });
        expect(cell(stats, 0, 1, 1)).toEqual({ mean: 3, frac: 100 });
    });

    it("only includes the given rows", () => {
        const stats = computeDotPlotStats({ ...base, rows: new Uint32Array([0, 2, 3, 4]) });
        expect(stats.panels[0].map((g) => g.count)).toEqual([2, 2]);
        expect(cell(stats, 0, 0, 0)).toEqual({ mean: 3, frac: 100 });
        expect(cell(stats, 0, 1, 0)).toEqual({ mean: 0.5, frac: 50 });
    });

    it("makes one panel per split value, in the given order, ignoring other values", () => {
        // split codes: 0 = Control, 1 = Disease, 2 = Treatment (not shown)
        const stats = computeDotPlotStats({
            ...base,
            split: { codes: new Uint8Array([0, 1, 0, 1, 2, 2]), panels: [1, 0] },
        });
        expect(stats.panels).toHaveLength(2);
        // Disease: data rows 1 and 3
        expect(stats.panels[0].map((g) => g.count)).toEqual([1, 1]);
        expect(cell(stats, 0, 0, 0)).toEqual({ mean: 0, frac: 0 });
        expect(cell(stats, 0, 1, 0)).toEqual({ mean: 1, frac: 100 });
        expect(cell(stats, 0, 1, 1)).toEqual({ mean: 3, frac: 100 });
        // Control: data rows 0 and 2
        expect(stats.panels[1].map((g) => g.count)).toEqual([2, 0]);
        expect(cell(stats, 1, 0, 0)).toEqual({ mean: 3, frac: 100 });
        expect(cell(stats, 1, 1, 0)).toEqual({ mean: 0, frac: 0 });
    });

    it("keeps every row category in every panel, including split values with no rows", () => {
        const stats = computeDotPlotStats({
            ...base,
            split: { codes: new Uint8Array([0, 0, 0, 0, 0, 0]), panels: [0, 5] },
        });
        expect(stats.panels[0].map((g) => g.count)).toEqual([3, 3]);
        expect(stats.panels[1]).toHaveLength(2);
        expect(stats.panels[1].map((g) => g.count)).toEqual([0, 0]);
    });

    it("gives one colour range across all panels", () => {
        const stats = computeDotPlotStats({
            rowCodes: new Uint8Array([0, 0, 0, 0]),
            rowCount: 1,
            fields: [new Float32Array([1, 1, 4, 4])],
            split: { codes: new Uint8Array([0, 0, 1, 1]), panels: [0, 1] },
        });
        expect(cell(stats, 0, 0, 0).mean).toBe(1);
        expect(cell(stats, 1, 0, 0).mean).toBe(4);
        expect(stats.meanRange).toEqual([1, 4]);
    });

    it("leaves empty groups out of the colour range", () => {
        const stats = computeDotPlotStats({
            rowCodes: new Uint8Array([0, 1]),
            rowCount: 3,
            fields: [new Float32Array([1, 2])],
        });
        expect(stats.panels[0][2]).toEqual({ count: 0, cells: [{ mean: 0, frac: 0 }] });
        expect(stats.meanRange).toEqual([1, 2]);
    });

    it("returns no colour range when every row is filtered out", () => {
        const stats = computeDotPlotStats({ ...base, rows: new Uint32Array(0) });
        expect(stats.meanRange).toBeNull();
        expect(stats.fracMax).toBe(0);
    });

    it("handles text16 codes above 255 for rows and splits", () => {
        const stats = computeDotPlotStats({
            rowCodes: new Uint16Array([300, 300, 5]),
            rowCount: 301,
            fields: [new Float32Array([1, 0, 2])],
            split: { codes: new Uint16Array([1000, 1000, 7]), panels: [1000] },
        });
        expect(stats.panels[0][300]).toEqual({ count: 2, cells: [{ mean: 0.5, frac: 50 }] });
        expect(stats.panels[0][5].count).toBe(0);
    });

    it("works with no fields", () => {
        const stats = computeDotPlotStats({ ...base, fields: [] });
        expect(stats.panels[0]).toEqual([
            { count: 3, cells: [] },
            { count: 3, cells: [] },
        ]);
        expect(stats.meanRange).toBeNull();
    });
});
