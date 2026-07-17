import { describe, expect, test } from "vitest";
import { formatHistogramBrushRange } from "@/react/components/HistogramWidget";

describe("formatHistogramBrushRange", () => {
    test("formats numeric ranges with compact g notation", () => {
        expect(formatHistogramBrushRange(0.5, 1.25)).toBe("0.5 - 1.25");
        expect(formatHistogramBrushRange(20150, 20320)).toBe("2.015e+4 - 2.032e+4");
    });

    test("formats date day ranges as YYYY-MM-DD", () => {
        // 2025-02-28 and 2025-08-23 as days since Unix epoch
        expect(formatHistogramBrushRange(20147, 20323, true)).toBe(
            "2025-02-28 - 2025-08-23",
        );
    });

    test("orders date range ends low-to-high", () => {
        expect(formatHistogramBrushRange(20323, 20147, true)).toBe(
            "2025-02-28 - 2025-08-23",
        );
    });
});
