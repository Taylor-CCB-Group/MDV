import { describe, expect, test } from "vitest";
import { formatDateDays } from "@/lib/dateFormat";

/**
 * Mirrors DataStore.getValue numeric+is_date branch without constructing a full DataStore
 * (Worker / config typing noise in unit tests).
 */
function getDateDisplayValue(raw: number, isDate: boolean): string | number {
    if (Number.isNaN(raw)) {
        return "missing";
    }
    if (isDate) {
        return formatDateDays(raw);
    }
    return raw;
}

describe("date column display / sort contract", () => {
    test("is_date columns display ISO dates while raw days stay numeric", () => {
        const raw = new Float32Array([0, 18262, Number.NaN]);
        expect(getDateDisplayValue(raw[0], true)).toBe("1970-01-01");
        expect(getDateDisplayValue(raw[1], true)).toBe("2020-01-01");
        expect(getDateDisplayValue(raw[2], true)).toBe("missing");
        expect(raw[0]).toBe(0);
        expect(raw[1]).toBe(18262);
    });

    test("numeric sort of day doubles is chronological", () => {
        const days = new Float32Array([18262, 0, 100]);
        const indices = new Uint32Array([0, 1, 2]);
        indices.sort((a, b) => {
            const comparison = days[a] < days[b] ? -1 : days[a] > days[b] ? 1 : 0;
            return comparison;
        });
        expect(Array.from(indices)).toEqual([1, 2, 0]);
    });
});
