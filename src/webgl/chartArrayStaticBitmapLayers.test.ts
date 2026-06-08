import { OrthographicViewport } from "@deck.gl/core";
import { describe, expect, test } from "vitest";
import {
    CHART_ARRAY_STATIC_VIEW_FLIP_Y,
    getFramebufferCompositeBounds,
    getOrthographicWorldBounds,
} from "./chartArrayStaticBitmapLayers";

describe("chartArrayStaticBitmapLayers", () => {
    test("framebuffer composite bounds match OrthographicViewport corners", () => {
        const viewState = {
            target: [100, 200, 0] as [number, number, number],
            zoom: 1,
        };
        const width = 400;
        const height = 300;
        const viewport = new OrthographicViewport({
            ...viewState,
            width,
            height,
            flipY: CHART_ARRAY_STATIC_VIEW_FLIP_Y,
            id: "test",
        });
        const topLeft = viewport.unproject([0, 0]);
        const bottomRight = viewport.unproject([width, height]);
        const corners = getFramebufferCompositeBounds(viewState, width, height);
        expect(corners[0][0]).toBeCloseTo(topLeft[0], 5);
        expect(corners[0][1]).toBeCloseTo(topLeft[1], 5);
        expect(corners[2][0]).toBeCloseTo(bottomRight[0], 5);
        expect(corners[2][1]).toBeCloseTo(bottomRight[1], 5);
    });

    test("flipY=true viewport inverts Y extent vs flipY=false formula", () => {
        const viewState = {
            target: [0, 0, 0] as [number, number, number],
            zoom: 0,
        };
        const width = 200;
        const height = 100;
        const [, bottom, , top] = getOrthographicWorldBounds(viewState, width, height);
        const corners = getFramebufferCompositeBounds(viewState, width, height);
        const compositeBottom = Math.min(corners[1][1], corners[2][1]);
        const compositeTop = Math.max(corners[0][1], corners[3][1]);
        expect(compositeBottom).not.toBeCloseTo(bottom, 1);
        expect(compositeTop).not.toBeCloseTo(top, 1);
    });
});
