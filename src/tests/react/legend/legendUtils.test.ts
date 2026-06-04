import { describe, expect, test } from "vitest";
import {
    formatLegendLabel,
    measureLegendLabelWidth,
} from "@/react/legend/shared/legendUtils";

describe("legendUtils", () => {
    test("does not truncate short labels", () => {
        const result = formatLegendLabel("B cell", 200);
        expect(result.truncated).toBe(false);
        expect(result.display).toBe("B cell");
        expect(result.full).toBe("B cell");
    });

    test("truncates long labels and preserves full text for tooltips", () => {
        const long =
            "CD8-positive cytotoxic T lymphocyte population cluster";
        const maxWidth = measureLegendLabelWidth("CD8-positive cy");
        const result = formatLegendLabel(long, maxWidth);
        expect(result.truncated).toBe(true);
        expect(result.display.endsWith("…")).toBe(true);
        expect(result.display.length).toBeLessThan(long.length);
        expect(result.full).toBe(long);
    });
});
