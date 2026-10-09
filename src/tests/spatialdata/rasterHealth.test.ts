import { describe, expect, it } from "vitest";
import {
    assessRasterLevels,
    loadRasterLevels,
    parseArrayMetadata,
    type RasterLevelInfo,
} from "@/react/spatialdata/raster_health";

const level = (path: string, shape: number[], chunks: number[], bytesPerElement = 2): RasterLevelInfo => ({
    path,
    shape,
    chunks,
    bytesPerElement,
});

describe("assessRasterLevels", () => {
    it("flags a single large level and full-width strip chunks (leap034 layout)", () => {
        const issues = assessRasterLevels([level("0", [37, 16333, 3788], [1, 4096, 3788], 4)]);
        expect(issues).toEqual([
            { kind: "no-pyramid", width: 3788, height: 16333 },
            { kind: "oversized-chunks", chunkWidth: 3788, chunkHeight: 4096, decodedBytes: 4096 * 3788 * 4 },
        ]);
    });

    it("accepts a pyramid that reaches a small coarsest level with tile-sized chunks", () => {
        const levels = [0, 1, 2, 3, 4].map((i) =>
            level(String(i), [4, Math.ceil(16333 / 2 ** i), Math.ceil(3788 / 2 ** i)], [1, 1024, 1024]),
        );
        expect(assessRasterLevels(levels)).toEqual([]);
    });

    it("flags a pyramid whose coarsest level is still large", () => {
        const issues = assessRasterLevels([
            level("0", [1, 40000, 40000], [1, 512, 512]),
            level("1", [1, 20000, 20000], [1, 512, 512]),
        ]);
        expect(issues).toEqual([{ kind: "coarsest-level-large", levels: 2, width: 20000, height: 20000 }]);
    });

    it("does not flag a small single-level image", () => {
        expect(assessRasterLevels([level("0", [3, 1000, 1000], [3, 1000, 1000], 1)])).toEqual([]);
    });
});

describe("parseArrayMetadata", () => {
    it("reads v3 zarr.json", () => {
        const meta = {
            shape: [4, 16333, 3788],
            data_type: "uint16",
            chunk_grid: { name: "regular", configuration: { chunk_shape: [1, 4096, 3788] } },
            codecs: [{ name: "bytes" }, { name: "zstd" }],
        };
        expect(parseArrayMetadata("0", meta)).toEqual(level("0", [4, 16333, 3788], [1, 4096, 3788], 2));
    });

    it("uses the inner chunk of a sharded v3 array", () => {
        const meta = {
            shape: [1, 8192, 8192],
            data_type: "float32",
            chunk_grid: { name: "regular", configuration: { chunk_shape: [1, 8192, 8192] } },
            codecs: [{ name: "sharding_indexed", configuration: { chunk_shape: [1, 512, 512] } }],
        };
        expect(parseArrayMetadata("0", meta)?.chunks).toEqual([1, 512, 512]);
    });

    it("reads v2 .zarray", () => {
        const meta = { shape: [3, 2048, 2048], chunks: [1, 512, 512], dtype: "<f4" };
        expect(parseArrayMetadata("1", meta)).toEqual(level("1", [3, 2048, 2048], [1, 512, 512], 4));
    });

    it("returns null for unrecognised metadata", () => {
        expect(parseArrayMetadata("0", undefined)).toBeNull();
        expect(parseArrayMetadata("0", { shape: [1], data_type: "string" })).toBeNull();
    });
});

describe("loadRasterLevels", () => {
    it("reads zarr.json per level, falling back to .zarray", async () => {
        const files: Record<string, unknown> = {
            "/0/zarr.json": {
                shape: [1, 100, 100],
                data_type: "uint8",
                chunk_grid: { configuration: { chunk_shape: [1, 100, 100] } },
            },
            "/1/.zarray": { shape: [1, 50, 50], chunks: [1, 50, 50], dtype: "|u1" },
        };
        const store = {
            get: async (key: string) =>
                key in files ? new TextEncoder().encode(JSON.stringify(files[key])) : undefined,
        };
        const levels = await loadRasterLevels({ scaleLevels: ["0", "1"], getStore: () => store });
        expect(levels.map((l) => [l.path, l.shape])).toEqual([
            ["0", [1, 100, 100]],
            ["1", [1, 50, 50]],
        ]);
    });
});
