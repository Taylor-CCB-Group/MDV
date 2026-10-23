import { computeDotPlotStats } from "@/react/components/dotplot/dotPlotStats";
import { describe, expect, it, vi } from "vitest";

type LegacyItem = { mean: number; frac: number };
type LegacyResult = { data: { count: number; values: LegacyItem[] }[]; frac_range: [number, number] };

/** Runs the legacy dot plot's `averages_simple` path in `catColWorker.js` on the calling thread. */
async function legacyAverages(
    rowCodes: Uint8Array,
    rowCount: number,
    fields: Float32Array[],
    filter: Uint8Array,
    threshold: number,
) {
    let result: LegacyResult | undefined;
    vi.stubGlobal("postMessage", (r: LegacyResult) => {
        result = r;
    });
    // @ts-expect-error a worker script, not a module: importing it installs its global onmessage
    await import("@/datastore/catColWorker.js");
    const handler = globalThis.onmessage as unknown as (e: { data: unknown }) => void;
    handler({
        data: [
            new Uint8Array(rowCodes.length).buffer, // the chart's own (local) filter: none
            filter.buffer,
            rowCodes.buffer,
            fields.map((f) => [f.buffer, "double"]),
            {
                method: "averages_simple",
                threshold,
                values: Array.from({ length: rowCount }, (_, i) => `cat${i}`),
                columns: fields.map((_, i) => `f${i}`),
                y_datatype: "text",
            },
        ],
    });
    vi.unstubAllGlobals();
    if (!result) throw new Error("legacy worker posted nothing");
    return result;
}

function randomData(seed: number, n: number, rowCount: number, fieldCount: number) {
    let s = seed;
    const rand = () => {
        s = (s * 1664525 + 1013904223) % 2 ** 32;
        return s / 2 ** 32;
    };
    const rowCodes = Uint8Array.from({ length: n }, () => Math.floor(rand() * rowCount));
    const fields = Array.from({ length: fieldCount }, () =>
        Float32Array.from({ length: n }, () => {
            const r = rand();
            if (r < 0.4) return 0;
            if (r < 0.45) return Number.NaN;
            return r * 5;
        }),
    );
    const filter = Uint8Array.from({ length: n }, () => (rand() < 0.2 ? 1 : 0));
    return { rowCodes, fields, filter };
}

describe("computeDotPlotStats matches the legacy dot plot", () => {
    it.each([0, 1.5])("with threshold %s", async (threshold) => {
        const rowCount = 6;
        const { rowCodes, fields, filter } = randomData(7, 2000, rowCount, 4);
        const legacy = await legacyAverages(rowCodes, rowCount, fields, filter, threshold);
        const rows = Uint32Array.from(filter.keys()).filter((i) => filter[i] === 0);
        const stats = computeDotPlotStats({ rowCodes, rowCount, fields, rows, threshold });
        for (let r = 0; r < rowCount; r++) {
            expect(stats.panels[0][r].count).toBe(legacy.data[r].count);
            for (let f = 0; f < fields.length; f++) {
                expect(stats.panels[0][r].cells[f].mean).toBeCloseTo(legacy.data[r].values[f].mean, 4);
                expect(stats.panels[0][r].cells[f].frac).toBeCloseTo(legacy.data[r].values[f].frac, 4);
            }
        }
        expect(stats.fracMax).toBeCloseTo(legacy.frac_range[1], 4);
    });
});
