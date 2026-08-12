import { describe, expect, test } from "vitest";
import {
    dateTickFormat,
    formatDateDays,
    isDateColumn,
    parseDateDays,
} from "@/lib/dateFormat";

describe("dateFormat", () => {
    test("formatDateDays formats UTC epoch days as YYYY-MM-DD", () => {
        expect(formatDateDays(0)).toBe("1970-01-01");
        // 2020-01-01 is 18262 days after 1970-01-01
        expect(formatDateDays(18262)).toBe("2020-01-01");
    });

    test("formatDateDays returns missing for non-finite", () => {
        expect(formatDateDays(Number.NaN)).toBe("missing");
        expect(formatDateDays(Number.POSITIVE_INFINITY)).toBe("missing");
    });

    test("parseDateDays parses YYYY-MM-DD to days", () => {
        expect(parseDateDays("1970-01-01")).toBe(0);
        expect(parseDateDays("2020-01-01")).toBe(18262);
        expect(parseDateDays("not-a-date")).toBeNull();
        expect(parseDateDays("")).toBeNull();
        expect(parseDateDays("18262")).toBeNull();
        expect(parseDateDays("2020-02-31")).toBeNull();
    });

    test("round-trip format and parse", () => {
        for (const days of [0, 1, 10000, 20000]) {
            expect(parseDateDays(formatDateDays(days))).toBe(days);
        }
    });

    test("isDateColumn and dateTickFormat", () => {
        expect(isDateColumn({ is_date: true, datatype: "double" })).toBe(true);
        expect(isDateColumn({ is_date: false, datatype: "double" })).toBe(false);
        expect(isDateColumn(null)).toBe(false);
        expect(dateTickFormat(0)).toBe("1970-01-01");
    });
});
