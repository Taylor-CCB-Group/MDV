import { describe, expect, test } from "vitest";
import { isCrosshairShown } from "./pointerCrosshair";

describe("pointerCrosshair", () => {
    test("marks the pointer unless the chart turns it off, and never in 3D", () => {
        expect(isCrosshairShown({})).toBe(true);
        expect(isCrosshairShown({ crosshair: null })).toBe(true);
        expect(isCrosshairShown({ crosshair: false })).toBe(false);
        expect(isCrosshairShown({ crosshair: true, dimension: "2d" })).toBe(true);
        expect(isCrosshairShown({ crosshair: true, dimension: "3d" })).toBe(false);
    });
});
