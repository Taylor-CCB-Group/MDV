import { describe, expect, test } from "vitest";
import {
    getChartArrayStaticCacheKey,
    getChartArrayStaticContentKey,
} from "./chartArrayRenderCache";

describe("chartArrayRenderCache", () => {
    test("buffer key rounds size and dpr only", () => {
        const keyA = getChartArrayStaticCacheKey({
            width: 1000.4,
            height: 700.2,
            devicePixelRatio: 2,
        });
        const keyB = getChartArrayStaticCacheKey({
            width: 1000.4,
            height: 700.2,
            devicePixelRatio: 2,
        });
        expect(keyA).toBe("1000x700@2");
        expect(keyA).toBe(keyB);
    });

    test("buffer key changes when size or dpr changes", () => {
        const base = getChartArrayStaticCacheKey({
            width: 1200,
            height: 800,
            devicePixelRatio: 1,
        });
        const widthChanged = getChartArrayStaticCacheKey({
            width: 1201,
            height: 800,
            devicePixelRatio: 1,
        });
        const dprChanged = getChartArrayStaticCacheKey({
            width: 1200,
            height: 800,
            devicePixelRatio: 2,
        });
        expect(base).not.toBe(widthChanged);
        expect(base).not.toBe(dprChanged);
    });

    test("content key is stable regardless of layer id order", () => {
        const viewState = { target: [10, 20, 0], zoom: 1.5 };
        const keyA = getChartArrayStaticContentKey(viewState, ["layer-b", "layer-a"]);
        const keyB = getChartArrayStaticContentKey(viewState, ["layer-a", "layer-b"]);
        expect(keyA).toBe(keyB);
    });

    test("content key changes when pan/zoom or layers change", () => {
        const base = getChartArrayStaticContentKey(
            { target: [0, 0, 0], zoom: 0 },
            ["scatter"],
        );
        const panned = getChartArrayStaticContentKey(
            { target: [1, 0, 0], zoom: 0 },
            ["scatter"],
        );
        const zoomed = getChartArrayStaticContentKey(
            { target: [0, 0, 0], zoom: 1 },
            ["scatter"],
        );
        const layersChanged = getChartArrayStaticContentKey(
            { target: [0, 0, 0], zoom: 0 },
            ["scatter-grey"],
        );
        expect(base).not.toBe(panned);
        expect(base).not.toBe(zoomed);
        expect(base).not.toBe(layersChanged);
    });
});
