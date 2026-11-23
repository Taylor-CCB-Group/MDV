import { describe, expect, test } from "vitest";
import { getViewportAtCanvasPoint } from "./deckViewports";

describe("getViewportAtCanvasPoint", () => {
    const viewports = [
        { id: "a", x: 0, y: 0, width: 100, height: 100 },
        { id: "b", x: 100, y: 0, width: 100, height: 100 },
    ];

    test("finds the viewport under the point", () => {
        expect(getViewportAtCanvasPoint(viewports, 10, 10)?.id).toBe("a");
        expect(getViewportAtCanvasPoint(viewports, 150, 10)?.id).toBe("b");
    });

    test("treats a viewport's right and bottom edges as outside it", () => {
        expect(getViewportAtCanvasPoint(viewports, 100, 10)?.id).toBe("b");
        expect(getViewportAtCanvasPoint(viewports, 10, 100)?.id).toBe("a");
    });

    test("prefers the viewport drawn last where they overlap", () => {
        const overlapping = [...viewports, { id: "c", x: 50, y: 0, width: 100, height: 100 }];
        expect(getViewportAtCanvasPoint(overlapping, 60, 10)?.id).toBe("c");
    });

    test("falls back to the first viewport off the canvas", () => {
        expect(getViewportAtCanvasPoint(viewports, 500, 500)?.id).toBe("a");
        expect(getViewportAtCanvasPoint([], 0, 0)).toBeUndefined();
    });
});
