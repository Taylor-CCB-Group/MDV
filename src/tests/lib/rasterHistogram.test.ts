import * as d3 from "d3";
import { describe, expect, it } from "vitest";
import { binRasterValues, computeRasterHistogram, sampleRasterValues } from "@/react/utils/rasterHistogram";

const sum = (values: number[]) => values.reduce((a, b) => a + b, 0);

describe("binRasterValues", () => {
    it("returns exactly bins counts and bins + 1 edges spanning the domain", () => {
        const { counts, edges } = binRasterValues(new Float64Array([0, 50, 100]), [0, 100], 200, "linear");
        expect(counts).toHaveLength(200);
        expect(edges).toHaveLength(201);
        expect(edges[0]).toBe(0);
        expect(edges[200]).toBe(100);
    });

    it("puts the domain max in the last bin and ignores out-of-domain and non-finite values", () => {
        const values = new Float64Array([0, 9.99, 10, -1, 11, Number.NaN, Number.POSITIVE_INFINITY]);
        const { counts } = binRasterValues(values, [0, 10], 10, "linear");
        expect(counts[0]).toBe(1);
        expect(counts[9]).toBe(2);
        expect(sum(counts)).toBe(3);
    });

    it("uses symlog-spaced edges matching d3.scaleSymlog in log mode", () => {
        const { edges } = binRasterValues(new Float64Array(), [0, 65535], 200, "log");
        const scale = d3.scaleSymlog().domain([0, 65535]).range([0, 1]);
        for (const i of [0, 1, 50, 100, 199, 200]) {
            expect(edges[i]).toBeCloseTo(scale.invert(i / 200), 6);
        }
    });

    it("counts each value into the log bin whose edges contain it", () => {
        const values = Float64Array.from({ length: 1000 }, (_, i) => i * 7);
        const { counts, edges } = binRasterValues(values, [0, 7000], 50, "log");
        expect(sum(counts)).toBe(1000);
        for (const value of values) {
            const bin = edges.findIndex((edge, i) => i < 50 && value >= edge && value < edges[i + 1]);
            expect(counts[bin === -1 ? 49 : bin]).toBeGreaterThan(0);
        }
    });

    it("widens a degenerate domain instead of dividing by zero", () => {
        const { counts, edges } = binRasterValues(new Float64Array([5, 5, 5]), [5, 5], 10, "linear");
        expect(counts[0]).toBe(3);
        expect(edges[10]).toBe(6);
    });
});

describe("sampleRasterValues", () => {
    it("copies small rasters verbatim", () => {
        const values = new Uint16Array([1, 2, 3, 4]);
        expect(Array.from(sampleRasterValues(values, 10))).toEqual([1, 2, 3, 4]);
    });

    it("caps large rasters at maxSamples, deterministically, without aliasing onto a column", () => {
        // A 1000-wide image where every value is its column index.
        const width = 1000;
        const values = Float32Array.from({ length: width * 1000 }, (_, i) => i % width);
        const sample = sampleRasterValues(values, 10_000);
        expect(sample).toHaveLength(10_000);
        expect(sampleRasterValues(values, 10_000)).toEqual(sample);
        // Stride is exactly 100, which divides the width: a fixed-offset sampler would see 10 columns.
        expect(new Set(sample).size).toBeGreaterThan(900);
    });
});

describe("computeRasterHistogram", () => {
    it("resolves auto x-scale from the values", () => {
        const skewed = Float64Array.from({ length: 1000 }, (_, i) => (i < 990 ? 1 : 60000));
        const result = computeRasterHistogram({ values: skewed, domain: [0, 65535], bins: 200, xScaleMode: "auto" });
        expect(result.xScale).toBe("log");
        const even = Float64Array.from({ length: 1000 }, (_, i) => i * 65);
        expect(computeRasterHistogram({ values: even, domain: [0, 65535], bins: 200, xScaleMode: "auto" }).xScale).toBe(
            "linear",
        );
    });
});
