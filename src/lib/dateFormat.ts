import type { DataColumn, DataType, NumberDataType } from "@/charts/charts";

const MS_PER_DAY = 86_400_000;

/** True when a column stores calendar dates as days since Unix epoch. */
export function isDateColumn(
    col: Pick<DataColumn<DataType>, "is_date" | "date_unit" | "datatype"> | null | undefined,
): boolean {
    return col?.is_date === true || col?.date_unit === "days";
}

/**
 * Format a day-number (days since Unix epoch, UTC) as `YYYY-MM-DD`.
 * Returns `"missing"` for non-finite values.
 */
export function formatDateDays(dayNumber: number): string {
    if (!Number.isFinite(dayNumber)) {
        return "missing";
    }
    const date = new Date(Math.round(dayNumber) * MS_PER_DAY);
    const y = date.getUTCFullYear();
    const m = String(date.getUTCMonth() + 1).padStart(2, "0");
    const d = String(date.getUTCDate()).padStart(2, "0");
    return `${y}-${m}-${d}`;
}

/**
 * Parse a strict ISO calendar date (`YYYY-MM-DD`) to days since Unix epoch
 * (UTC midnight). Returns `null` if invalid. Bare numbers and loose Date.parse
 * forms are rejected so values like `"18262"` are not treated as years.
 */
export function parseDateDays(iso: string): number | null {
    const trimmed = iso.trim();
    if (!trimmed) {
        return null;
    }
    const m = /^(\d{4})-(\d{2})-(\d{2})$/.exec(trimmed);
    if (!m) {
        return null;
    }
    const year = Number(m[1]);
    const month = Number(m[2]);
    const day = Number(m[3]);
    const ms = Date.UTC(year, month - 1, day);
    // Reject impossible calendar dates (e.g. 2020-02-31 → Mar 2).
    const date = new Date(ms);
    if (
        date.getUTCFullYear() !== year ||
        date.getUTCMonth() + 1 !== month ||
        date.getUTCDate() !== day
    ) {
        return null;
    }
    return Math.round(ms / MS_PER_DAY);
}

/** Tick formatter for continuous axes bound to an `is_date` column. */
export function dateTickFormat(value: number): string {
    return formatDateDays(value);
}

export type DateAwareNumberColumn = DataColumn<NumberDataType>;
