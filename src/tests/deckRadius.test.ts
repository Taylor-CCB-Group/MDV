import { expect, test } from "vitest";
import { deckRadiusForPixelRadius, radiusReferenceRange3D } from "@/react/axisScaling";

const total = (s: { radius: number; course_radius: number } | undefined) => (s ? s.radius * s.course_radius : 0);

test("a point size beyond the slider's range is kept, not capped", () => {
    // a classic 2.5px radius at a view where a data unit is tiny on screen
    const sizing = deckRadiusForPixelRadius(2.5, 0.15, 1);
    expect(total(sizing)).toBeCloseTo((2.5 * 10000) / 0.15, -2);
    expect(sizing?.course_radius).toBe(100);
});

test("sizes that fit the slider use the coarsest setting they need", () => {
    expect(deckRadiusForPixelRadius(2.5, 1000, 100)).toEqual({ radius: 0.3, course_radius: 1 });
    expect(deckRadiusForPixelRadius(2.5, 100, 100)).toEqual({ radius: 2.5, course_radius: 1 });
    expect(deckRadiusForPixelRadius(2.5, 10, 100)).toEqual({ radius: 2.5, course_radius: 10 });
});

test("a 3D chart's points are sized from its largest range", () => {
    expect(radiusReferenceRange3D([0, 0, 16])).toBe(16);
    expect(radiusReferenceRange3D([40, 20, 16])).toBe(40);
});
