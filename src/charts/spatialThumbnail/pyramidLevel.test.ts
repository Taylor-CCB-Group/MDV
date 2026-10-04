import { describe, expect, it } from "vitest";
import { compositePlanes, pickPyramidLevel, planeDomain, scalePlane } from "./pyramidLevel";

describe("pickPyramidLevel", () => {
    const levels = [
        { width: 1024, height: 1024 },
        { width: 256, height: 256 },
        { width: 64, height: 64 },
    ];

    it("picks the finest level that fits", () => {
        expect(pickPyramidLevel(levels, 128)).toBe(2);
        expect(pickPyramidLevel(levels, 256)).toBe(1);
        expect(pickPyramidLevel(levels, 64)).toBe(2);
    });

    it("uses the coarsest level when every level is larger", () => {
        expect(pickPyramidLevel(levels, 32)).toBe(2);
        expect(pickPyramidLevel([{ width: 400, height: 200 }], 128)).toBe(0);
    });

    it("prefers an earlier level when several fit", () => {
        expect(pickPyramidLevel(levels, 2048)).toBe(0);
    });
});

describe("scalePlane", () => {
    it("samples the source into a square", () => {
        const src = new Float32Array([1, 2, 3, 4]);
        const scaled = scalePlane(src, 2, 2, 2);
        expect(Array.from(scaled)).toEqual([1, 2, 3, 4]);
    });
});

describe("compositePlanes", () => {
    it("colors a visible channel and skips a hidden one", () => {
        const plane = new Float32Array([10]);
        const rgba = compositePlanes(1, [
            {
                plane,
                color: [0, 255, 0],
                visible: false,
                contrastLimits: [0, 10],
                brightness: 0.5,
                contrast: 0.5,
            },
        ]);
        expect(rgba[1]).toBe(0);
        const shown = compositePlanes(1, [
            {
                plane,
                color: [0, 255, 0],
                visible: true,
                contrastLimits: [0, 10],
                brightness: 0.5,
                contrast: 0.5,
            },
        ]);
        expect(shown[1]).toBeGreaterThan(0);
        expect(shown[0]).toBe(0);
    });
});

describe("planeDomain", () => {
    it("returns the finite min and max", () => {
        expect(planeDomain(new Float32Array([2, Number.NaN, 8, 3]))).toEqual([2, 8]);
    });
});
