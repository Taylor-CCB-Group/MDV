import type { SpatialData } from "@spatialdata/core";
import { describe, expect, it, vi } from "vitest";
import { createSpatialDataReader } from "@/react/spatialdata/spatial_data_reader";

const fakeSpatialData = (url: string) => ({ url }) as unknown as SpatialData;

describe("createSpatialDataReader", () => {
    it("reads a store once however many trees ask for it", async () => {
        const load = vi.fn(async (url: string) => fakeSpatialData(url));
        const reader = createSpatialDataReader(load);

        const chart = reader.read("a");
        await chart.promise;
        const dialog = reader.read("a");

        expect(load).toHaveBeenCalledTimes(1);
        expect(dialog).toBe(chart);
    });

    it("has the result ready for a tree that mounts after the read settled", async () => {
        const reader = createSpatialDataReader(async (url) => fakeSpatialData(url));
        const first = reader.read("a");
        expect(first.state).toEqual({ spatialData: null, loading: true, error: null });

        const spatialData = await first.promise;

        expect(reader.read("a").state).toEqual({ spatialData, loading: false, error: null });
    });

    it("keeps one read per URL, so trees disagreeing about the URL don't evict each other", () => {
        const load = vi.fn(async (url: string) => fakeSpatialData(url));
        const reader = createSpatialDataReader(load);

        const a = reader.read("a");
        const b = reader.read("b");

        expect(b).not.toBe(a);
        expect(reader.read("a")).toBe(a);
        expect(reader.read("b")).toBe(b);
        expect(load).toHaveBeenCalledTimes(2);
    });

    it("keeps a failure rather than retrying on every render", async () => {
        const load = vi.fn(async () => {
            throw "404";
        });
        const reader = createSpatialDataReader(load);

        await reader.read("a").promise.catch(() => {});
        const { state } = reader.read("a");

        expect(load).toHaveBeenCalledTimes(1);
        expect(state.loading).toBe(false);
        expect(state.error).toBeInstanceOf(Error);
        expect(state.error?.message).toBe("404");
    });

    it("reads a failed store again on retry and tells subscribers", async () => {
        const load = vi
            .fn()
            .mockRejectedValueOnce(new Error("timeout"))
            .mockImplementation(async (url: string) => fakeSpatialData(url));
        const reader = createSpatialDataReader(load);
        const listener = vi.fn();
        reader.subscribe(listener);

        await reader.read("a").promise.catch(() => {});
        listener.mockClear();
        reader.retry();

        expect(listener).toHaveBeenCalledTimes(1);
        expect(reader.read("a").state.loading).toBe(true);
        const spatialData = await reader.read("a").promise;
        expect(reader.read("a").state).toEqual({ spatialData, loading: false, error: null });
        expect(load).toHaveBeenCalledTimes(2);
    });

    it("ignores retry while the read is pending or succeeded", async () => {
        const load = vi.fn(async (url: string) => fakeSpatialData(url));
        const reader = createSpatialDataReader(load);

        reader.read("a");
        reader.retry();
        await reader.read("a").promise;
        reader.retry();

        expect(load).toHaveBeenCalledTimes(1);
    });
});
