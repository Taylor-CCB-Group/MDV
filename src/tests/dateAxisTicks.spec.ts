import { describe, expect, test } from "vitest";
import { dateTickFormat, formatDateDays } from "@/lib/dateFormat";

describe("date axis tick labels", () => {
    test("formats typical qc_date day ticks as ISO dates", () => {
        // Values visible on the MICRON-scanr-sim-psf scatter axis before formatting.
        expect(dateTickFormat(20100)).toBe("2025-01-12");
        expect(dateTickFormat(20150)).toBe("2025-03-03");
        expect(dateTickFormat(20454)).toBe(formatDateDays(20454));
        expect(dateTickFormat(20454)).toBe("2026-01-01");
    });

    test("d3 NumberValue-like objects coerce correctly", () => {
        const asNumberValue = { valueOf: () => 20089 };
        expect(dateTickFormat(Number(asNumberValue))).toBe("2025-01-01");
    });
});
