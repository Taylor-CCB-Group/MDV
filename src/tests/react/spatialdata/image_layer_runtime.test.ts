import { describe, expect, test } from "vitest";

import {
    channelOptionsForSelections,
    createStableChannelId,
    patchVivToneArrays,
    projectRuntimeCacheToArrays,
    pruneRuntimeCache,
    readVivToneArrays,
    syncViewerChannelArraysFromStore,
} from "@/react/spatialdata/image_layer_runtime";

describe("image_layer_runtime", () => {
    test("syncViewerChannelArraysFromStore labels channels by selection after middle remove", () => {
        const synced = syncViewerChannelArraysFromStore(
            {
                ids: ["a", "c"],
                colors: [[1, 0, 0], [0, 0, 1]],
                contrastLimits: [
                    [0, 100],
                    [50, 200],
                ],
                channelsVisible: [true, true],
                selections: [
                    { z: 0, c: 0, t: 0 },
                    { z: 0, c: 2, t: 0 },
                ],
            },
            {
                channelOptions: ["DAPI", "GFP", "Cy5"],
                isChannelLoading: [false, true, false],
                pixelValues: [1, 2, 3],
            },
            ["DAPI", "GFP", "Cy5"],
        );
        expect(synced.channelOptions).toEqual(["DAPI", "Cy5"]);
        expect(synced.isChannelLoading).toEqual([false, false]);
        expect(synced.pixelValues.every((value) => Number.isNaN(value))).toBe(true);
    });

    test("projectRuntimeCacheToArrays preserves trailing channel stats by id", () => {
        const cache = new Map([
            ["a", { domains: [0, 10] as [number, number], raster: { width: 1, height: 1, data: new Float32Array([1]) } }],
            ["c", { domains: [40, 50] as [number, number], raster: { width: 1, height: 1, data: new Float32Array([3]) } }],
        ]);
        const projected = projectRuntimeCacheToArrays(
            ["a", "c"],
            cache,
            [
                [0, 100],
                [50, 200],
            ],
        );
        expect(projected.domains).toEqual([
            [0, 10],
            [40, 50],
        ]);
        expect(projected.raster[1]?.data[0]).toBe(3);
    });

    test("pruneRuntimeCache drops removed ids", () => {
        const cache = new Map([
            ["a", { domains: [0, 1] as [number, number], raster: { width: 0, height: 0, data: new Float32Array() } }],
            ["b", { domains: [0, 1] as [number, number], raster: { width: 0, height: 0, data: new Float32Array() } }],
        ]);
        pruneRuntimeCache(cache, ["a"]);
        expect(cache.has("b")).toBe(false);
        expect(cache.has("a")).toBe(true);
    });

    test("readVivToneArrays and patchVivToneArrays round-trip", () => {
        const patched = patchVivToneArrays(undefined, 2, 1, { brightness: 0.2, contrast: 0.8 });
        expect(readVivToneArrays(patched, 2)).toEqual({
            brightness: [0.5, 0.2],
            contrast: [0.5, 0.8],
        });
    });

    test("createStableChannelId includes layer prefix", () => {
        expect(createStableChannelId("spatialdata-image-a")).toMatch(/^spatialdata-image-a-ch-.+/);
    });

    test("channelOptionsForSelections uses selection.c", () => {
        expect(
            channelOptionsForSelections(
                [
                    { z: 0, c: 0, t: 0 },
                    { z: 0, c: 2, t: 0 },
                ],
                ["DAPI", "GFP", "Cy5"],
            ),
        ).toEqual(["DAPI", "Cy5"]);
    });
});
