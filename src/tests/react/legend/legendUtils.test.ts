import { describe, expect, test } from "vitest";
import {
    formatLegendLabel,
    getContinuousLegendContainerHeight,
    getContinuousLegendLayout,
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

    test("continuous layout leaves room for axis and labels", () => {
        const layout = getContinuousLegendLayout(180, true);

        expect(layout.barX).toBeGreaterThan(0);
        expect(layout.axisWidth).toBeGreaterThan(0);
        expect(layout.axisWidth).toBeLessThan(layout.layoutWidth);
        expect(layout.labelMaxWidth).toBeGreaterThan(0);
        expect(layout.tickCount).toBeGreaterThanOrEqual(2);
    });

    test("continuous legend height grows for longer tick labels", () => {
        const options = { width: 180, hasLabel: true };
        expect(
            getContinuousLegendContainerHeight([0, 10], options),
        ).toBeLessThan(
            getContinuousLegendContainerHeight(
                [-0.123456789, 123456789],
                options,
            ),
        );
    });

    test("continuous legend height fits rotated scientific notation ticks", () => {
        const layout = getContinuousLegendLayout(180, true);
        const height = getContinuousLegendContainerHeight([0, 1e10], {
            width: 180,
            hasLabel: true,
        });

        expect(height).toBeGreaterThan(
            layout.barY + layout.barHeight + 30,
        );
        expect(height).toBeGreaterThan(
            getContinuousLegendContainerHeight([0, 10], {
                width: 180,
                hasLabel: true,
            }),
        );
    });
});
