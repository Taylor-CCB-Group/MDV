import type { SpatialData } from "@spatialdata/core";
import { createSpatialDataReader } from "@/react/spatialdata/spatial_data_reader";
import { describe, expect, it, vi } from "vitest";

const fakeSpatialData = (url: string) => ({ url }) as unknown as SpatialData;

describe("createSpatialDataReader", () => {
    it("reads a store once however many trees ask for it", async () => {
        const read = vi.fn(async (url: string) => fakeSpatialData(url));
        const readSpatialData = createSpatialDataReader(read);

        const chart = readSpatialData("a");
        await chart.promise;
        const dialog = readSpatialData("a");

        expect(read).toHaveBeenCalledTimes(1);
        expect(dialog).toBe(chart);
    });

    it("has the result ready for a tree that mounts after the read settled", async () => {
        const readSpatialData = createSpatialDataReader(async (url) => fakeSpatialData(url));
        const first = readSpatialData("a");
        expect(first.state).toEqual({ spatialData: null, loading: true, error: null });

        const spatialData = await first.promise;

        expect(readSpatialData("a").state).toEqual({ spatialData, loading: false, error: null });
    });

    it("reads again when the URL changes", () => {
        const read = vi.fn(async (url: string) => fakeSpatialData(url));
        const readSpatialData = createSpatialDataReader(read);

        const a = readSpatialData("a");
        const b = readSpatialData("b");

        expect(read).toHaveBeenCalledTimes(2);
        expect(b).not.toBe(a);
    });

    it("keeps a failure rather than retrying on every render", async () => {
        const read = vi.fn(async () => {
            throw "404";
        });
        const readSpatialData = createSpatialDataReader(read);

        await readSpatialData("a").promise.catch(() => {});
        const { state } = readSpatialData("a");

        expect(read).toHaveBeenCalledTimes(1);
        expect(state.loading).toBe(false);
        expect(state.error).toBeInstanceOf(Error);
        expect(state.error?.message).toBe("404");
    });
});
