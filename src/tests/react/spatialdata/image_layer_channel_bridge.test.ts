import { describe, expect, test } from "vitest";

import {
    buildChannelDefaultsFromConfig,
    channelStateKey,
    createStableChannelId,
    mergeHydrateChannelsState,
    persistedChannelCount,
    serializeChannelsFromStore,
    syncViewerChannelArraysFromStore,
} from "@/react/spatialdata/image_layer_channel_bridge";

describe("image_layer_channel_bridge", () => {
    test("persistedChannelCount ignores loaderDefaults when channels exist", () => {
        const count = persistedChannelCount(
            {
                channelIds: ["a", "b"],
                colors: [[255, 0, 0], [0, 255, 0]],
            },
            {
                colors: [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
            },
            ["Ch0", "Ch1", "Ch2"],
        );
        expect(count).toBe(2);
    });

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

    test("mergeHydrateChannelsState preserves domains by stable id", () => {
        const defaults = buildChannelDefaultsFromConfig({
            layerId: "layer",
            channels: {
                channelIds: ["a", "c"],
                contrastLimits: [
                    [0, 100],
                    [50, 200],
                ],
                selections: [{ c: 0 }, { c: 2 }],
            },
            loaderDefaults: undefined,
            channelNames: ["DAPI", "GFP", "Cy5"],
        });
        const merged = mergeHydrateChannelsState(defaults, {
            ids: ["a", "c"],
            domains: [
                [0, 10],
                [40, 50],
            ],
            raster: [
                { width: 1, height: 1, data: new Float32Array([1]) },
                { width: 1, height: 1, data: new Float32Array([3]) },
            ],
        });
        expect(merged.domains).toEqual([
            [0, 10],
            [40, 50],
        ]);
        expect(merged.raster[1]?.data[0]).toBe(3);
    });

    test("middle remove keeps trailing channel state by stable id", () => {
        const channels = {
            channelIds: ["layer-ch-aaa", "layer-ch-ccc"],
            colors: [
                [255, 0, 0],
                [0, 0, 255],
            ] as [number, number, number][],
            contrastLimits: [
                [0, 100],
                [50, 200],
            ] as [number, number][],
            channelsVisible: [true, false],
            selections: [{ c: 0 }, { c: 2 }],
        };

        const defaults = buildChannelDefaultsFromConfig({
            layerId: "layer",
            channels,
            loaderDefaults: {
                colors: [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                selections: [{ c: 0 }, { c: 1 }, { c: 2 }],
            },
            channelNames: ["DAPI", "GFP", "Cy5"],
        });

        expect(defaults.ids).toEqual(["layer-ch-aaa", "layer-ch-ccc"]);
        expect(defaults.colors[1]).toEqual([0, 0, 255]);
        expect(defaults.contrastLimits[1]).toEqual([50, 200]);
        expect(defaults.selections[1]?.c).toBe(2);
        expect(defaults.channelsVisible[1]).toBe(false);
    });

    test("serializeChannelsFromStore preserves ids without reindexing", () => {
        const serialized = serializeChannelsFromStore(
            {
                ids: ["layer-ch-one", "layer-ch-two"],
                colors: [[1, 2, 3], [4, 5, 6]],
                contrastLimits: [
                    [0, 10],
                    [20, 30],
                ],
                channelsVisible: [true, true],
                selections: [
                    { z: 0, c: 0, t: 0 },
                    { z: 0, c: 2, t: 0 },
                ],
            },
            {},
        );
        expect(serialized.channelIds).toEqual(["layer-ch-one", "layer-ch-two"]);
        expect(serialized.selections?.[1]?.c).toBe(2);
    });

    test("channelStateKey changes when middle channel removed", () => {
        const before = channelStateKey({
            channelIds: ["a", "b", "c"],
            colors: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            channelsVisible: [true, true, true],
            selections: [{ c: 0 }, { c: 1 }, { c: 2 }],
        });
        const after = channelStateKey({
            channelIds: ["a", "c"],
            colors: [[1, 0, 0], [0, 0, 1]],
            channelsVisible: [true, true],
            selections: [{ c: 0 }, { c: 2 }],
        });
        expect(before).not.toBe(after);
    });

    test("createStableChannelId includes layer prefix", () => {
        expect(createStableChannelId("spatialdata-image-a")).toMatch(
            /^spatialdata-image-a-ch-.+/,
        );
    });
});
