import { scaleLinear } from "d3-scale";
import { describe, expect, it, vi } from "vitest";
import { computeMultiLineDensities } from "@/react/components/multiline/multiLineStats";

type LegacyCurve = number[] & { id: number; max: number };

/** Runs the legacy multi-line chart's `kernal_density` path in `boxPlotWorker.js` on the calling thread. */
async function legacyDensities(
    categoryCodes: Uint8Array,
    categoryCount: number,
    values: Float32Array,
    filter: Uint8Array,
    config: { ticks: number[]; bandwidth: number; scaletrim?: [number, number] },
) {
    let result: LegacyCurve[] | undefined;
    vi.stubGlobal("postMessage", (r: LegacyCurve[]) => {
        result = r;
    });
    // @ts-expect-error a worker script, not a module: importing it installs its global onmessage
    await import("@/datastore/boxPlotWorker.js");
    const handler = globalThis.onmessage as unknown as (e: { data: unknown }) => void;
    handler({
        data: [
            new Uint8Array(values.length).buffer, // the chart's own (local) filter: none
            filter.buffer,
            categoryCodes.buffer,
            [values.buffer, "double"],
            {
                type: "kernal_density",
                values: Array.from({ length: categoryCount }, (_, i) => `cat${i}`),
                xPos: new Float32Array(values.length).buffer,
                ...config,
            },
        ],
    });
    vi.unstubAllGlobals();
    if (!result) throw new Error("legacy worker posted nothing");
    return result;
}

/** A small deterministic dataset: three categories with different spreads, some NaNs and filtered rows. */
function makeData(n = 500) {
    let seed = 1;
    const random = () => {
        seed = (seed * 16807) % 2147483647;
        return (seed - 1) / 2147483646;
    };
    const categoryCodes = new Uint8Array(n);
    const values = new Float32Array(n);
    const filter = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
        const c = i % 3;
        categoryCodes[i] = c;
        values[i] = c * 2 + (random() - 0.5) * (c + 1) * 3;
        if (i % 17 === 0) values[i] = Number.NaN;
        if (i % 11 === 0) filter[i] = 1;
    }
    const rows = Uint32Array.from({ length: n }, (_, i) => i).filter((i) => filter[i] === 0);
    return { categoryCodes, values, filter, rows };
}

describe("computeMultiLineDensities", () => {
    it("matches the legacy worker exactly without a trim", async () => {
        const { categoryCodes, values, filter, rows } = makeData();
        const ticks = scaleLinear().domain([-6, 10]).ticks(40);
        const bandwidth = 0.7;
        const legacy = await legacyDensities(categoryCodes, 4, values, filter, { ticks, bandwidth });
        const { curves } = computeMultiLineDensities({
            categoryCodes,
            categoryCount: 4,
            values,
            rows,
            ticks,
            bandwidth,
        });
        // category 3 has no rows: neither returns a curve for it
        expect(curves.map((c) => c.category)).toEqual(legacy.map((c) => c.id));
        curves.forEach((curve, i) => {
            expect(Array.from(curve.density)).toEqual(Array.from(legacy[i]));
            expect(curve.max).toBe(legacy[i].max);
        });
    });

    it("matches the legacy trimmed shape when clamping in data units", async () => {
        const { categoryCodes, values, filter, rows } = makeData();
        const [q0, q1] = [-1, 5];
        const span = q1 - q0;
        // legacy works on values normalised to [0, 1] of the trimmed range
        const legacyTicks = scaleLinear().domain([0, 1]).ticks(40);
        const legacy = await legacyDensities(categoryCodes, 3, values, filter, {
            ticks: legacyTicks,
            bandwidth: 0.1,
            scaletrim: [q0, q1],
        });
        const { curves } = computeMultiLineDensities({
            categoryCodes,
            categoryCount: 3,
            values,
            rows,
            ticks: legacyTicks.map((t) => q0 + t * span),
            bandwidth: 0.1 * span,
            clamp: [q0, q1],
        });
        expect(curves.map((c) => c.category)).toEqual(legacy.map((c) => c.id));
        curves.forEach((curve, i) => {
            // same curve, as a density per data unit rather than per normalised unit
            Array.from(curve.density).forEach((d, j) => {
                expect(d * span).toBeCloseTo(legacy[i][j], 10);
            });
        });
    });

    it("leaves out categories whose values are all NaN", () => {
        const { curves } = computeMultiLineDensities({
            categoryCodes: new Uint8Array([0, 0, 1, 1]),
            categoryCount: 2,
            values: new Float32Array([1, 2, Number.NaN, Number.NaN]),
            ticks: [0, 1, 2, 3],
            bandwidth: 1,
        });
        expect(curves.map((c) => c.category)).toEqual([0]);
    });

    it("reads text16 category codes", () => {
        const { curves } = computeMultiLineDensities({
            categoryCodes: new Uint16Array([300, 300, 2]),
            categoryCount: 301,
            values: new Float32Array([1, 1, 1]),
            ticks: [0, 1, 2],
            bandwidth: 1,
        });
        expect(curves.map((c) => c.category)).toEqual([2, 300]);
    });
});
