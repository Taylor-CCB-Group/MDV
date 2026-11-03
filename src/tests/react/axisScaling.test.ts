import { deckRadiusForClassicRadius, fitOrthographicView, radiusReferenceRange, serialisableZoom } from "@/react/axisScaling";
import { describe, expect, test } from "vitest";

describe("fitOrthographicView", () => {
    const wide = { minX: 0, maxX: 1, minY: 0, maxY: 10000 };

    test("equal uses one zoom, set by the axis that needs more room", () => {
        const view = fitOrthographicView(wide, 500, 400, "equal");
        expect(view.zoomX).toBeUndefined();
        expect(view.zoomY).toBeUndefined();
        expect(2 ** view.zoom * 10000).toBeLessThanOrEqual(400);
        // x is a sliver: far less than the chart width
        expect(2 ** view.zoom * 1).toBeLessThan(5);
    });

    test("fit zooms each axis to fill the chart", () => {
        const view = fitOrthographicView(wide, 500, 400, "fit");
        expect(view.zoomX).toBeDefined();
        expect(view.zoomY).toBeDefined();
        // both ranges take up (nearly) the same fraction of their axis
        const fractionX = (2 ** (view.zoomX ?? 0) * 1) / 500;
        const fractionY = (2 ** (view.zoomY ?? 0) * 10000) / 400;
        expect(fractionX).toBeCloseTo(fractionY);
        expect(fractionX).toBeGreaterThan(0.5);
        expect(view.zoom).toBe(Math.min(view.zoomX ?? 0, view.zoomY ?? 0));
    });

    test("centres on the data", () => {
        expect(fitOrthographicView({ minX: 10, maxX: 20, minY: -4, maxY: 0 }, 300, 300, "fit").target).toEqual([15, -2, 0]);
    });

    test("equal and fit agree when the ranges match the chart shape", () => {
        const bounds = { minX: 0, maxX: 100, minY: 0, maxY: 100 };
        const equal = fitOrthographicView(bounds, 400, 400, "equal");
        const fit = fitOrthographicView(bounds, 400, 400, "fit");
        expect(fit.zoomX).toBeCloseTo(equal.zoom);
        expect(fit.zoomY).toBeCloseTo(equal.zoom);
    });

    test("a single point gives a finite view", () => {
        const view = fitOrthographicView({ minX: 5, maxX: 5, minY: 5, maxY: 5 }, 300, 300, "fit");
        expect(Number.isFinite(view.zoom)).toBe(true);
        expect(Number.isFinite(view.zoomX)).toBe(true);
    });
});

describe("radiusReferenceRange", () => {
    test("equal is the larger range", () => {
        expect(radiusReferenceRange(1, 10000, { width: 500, height: 400 }, "equal")).toBe(10000);
    });

    test("fit makes points the same size on screen however different the ranges are", () => {
        const chart = { width: 500, height: 400 };
        const pixelRadius = (xRange: number, yRange: number) => {
            const bounds = { minX: 0, maxX: xRange, minY: 0, maxY: yRange };
            const { zoom } = fitOrthographicView(bounds, chart.width, chart.height, "fit");
            // deck sizes points by the viewport's (smaller) zoom
            return radiusReferenceRange(xRange, yRange, chart, "fit") * 2 ** zoom;
        };
        expect(pixelRadius(1, 10000)).toBeCloseTo(pixelRadius(500, 400));
        expect(pixelRadius(10000, 1)).toBeCloseTo(pixelRadius(500, 400));
    });

    test("fit with a square chart and square data matches equal", () => {
        expect(radiusReferenceRange(50, 50, { width: 300, height: 300 }, "fit")).toBeCloseTo(50);
    });

    test("falls back to equal for a chart with no size", () => {
        expect(radiusReferenceRange(3, 7, { width: 0, height: 0 }, "fit")).toBe(7);
    });
});

describe("serialisableZoom", () => {
    test("reduces a zoom pair to a number", () => {
        expect(serialisableZoom({ zoom: [2, 5], zoomX: 2, zoomY: 5 }, "equal")).toEqual({ zoom: 2 });
    });

    test("keeps both axes for fit", () => {
        expect(serialisableZoom({ zoom: [2, 5], zoomX: 2, zoomY: 5 }, "fit")).toEqual({ zoom: 2, zoomX: 2, zoomY: 5 });
    });

    test("fit without axis zooms is a plain zoom", () => {
        expect(serialisableZoom({ zoom: 3 }, "fit")).toEqual({ zoom: 3 });
    });
});

describe("deckRadiusForClassicRadius", () => {
    test("draws a point as large on screen as the classic chart did", () => {
        const size: [number, number] = [493, 462];
        const sizing = deckRadiusForClassicRadius(2, size);
        expect(sizing).toBeDefined();
        const { radius, course_radius } = sizing!;
        // the fitted view of data as wide as it is tall, in a square the chart's shorter side
        const view = fitOrthographicView({ minX: 0, maxX: 10, minY: 0, maxY: 10 }, size[1], size[1], "equal");
        const radiusPx = ((radius * course_radius) / 10000) * 10 * 2 ** view.zoom;
        expect(radiusPx * 2).toBeCloseTo(2, 1);
    });

    test("keeps the slider's range by moving large values into course_radius", () => {
        expect(deckRadiusForClassicRadius(2, [493, 462])).toEqual({ radius: 3.3, course_radius: 10 });
        expect(deckRadiusForClassicRadius(2, [8000, 8000])?.course_radius).toBe(1);
    });

    test("leaves the setting alone for a radius it cannot map", () => {
        expect(deckRadiusForClassicRadius(0, [300, 300])).toBeUndefined();
        expect(deckRadiusForClassicRadius(Number.NaN, [300, 300])).toBeUndefined();
    });
});

describe("fitOrthographicView with a constant column", () => {
    test("sizes the minimum span of each axis from its own bounds", () => {
        const view = fitOrthographicView({ minX: 1e6, maxX: 1e6, minY: 1, maxY: 1 }, 400, 400, "fit");
        // y is a constant near 1: its span is about 1e-6, not the 1 that x's bounds would give
        expect(view.zoomY).toBeGreaterThan((view.zoomX ?? 0) + 10);
    });
});
