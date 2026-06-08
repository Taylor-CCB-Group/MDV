import { describe, expect, test } from "vitest";
import {
    indicesFromIntersectionEntries,
    roundRect,
    snapRectToDevicePixels,
} from "./chartArrayUtils";

describe("chartArrayUtils", () => {
    test("snapRectToDevicePixels aligns bounds to the device pixel grid", () => {
        expect(snapRectToDevicePixels({ x: 10.3, y: 0, width: 50.4, height: 40.2 }, 2)).toEqual({
            x: 10.5,
            y: 0,
            width: 50,
            height: 40,
        });
    });

    test("roundRect snaps edges to integer pixels", () => {
        expect(roundRect({ x: 12.4, y: 8.6, width: 100.2, height: 99.7 })).toEqual({
            x: 12,
            y: 9,
            width: 101,
            height: 99,
        });
    });

    test("indicesFromIntersectionEntries returns visible sorted indices", () => {
        const target = document.createElement("div");
        target.dataset.chartArrayIndex = "2";
        const entries = [
            { target, isIntersecting: true },
        ] as unknown as IntersectionObserverEntry[];

        expect(indicesFromIntersectionEntries(entries, 5)).toEqual([2]);
    });

    test("indicesFromIntersectionEntries ignores out-of-range indices", () => {
        const target = document.createElement("div");
        target.dataset.chartArrayIndex = "9";
        const entries = [
            { target, isIntersecting: true },
        ] as unknown as IntersectionObserverEntry[];

        expect(indicesFromIntersectionEntries(entries, 3)).toEqual([]);
    });
});
