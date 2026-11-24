import * as d3 from "d3";
import { formatDateDays } from "@/lib/dateFormat";
import { useTheme } from "../../hooks";
import type { BrushXStyle } from "./useBrushX";

const formatNumericBrushValue = d3.format(".4~g");

/** Brush chip text for a numeric or calendar-date range. */
export function formatBrushRange(start: number, end: number, isDate = false): string {
    const lo = start <= end ? start : end;
    const hi = start <= end ? end : start;
    if (isDate) {
        return `${formatDateDays(lo)} - ${formatDateDays(hi)}`;
    }
    return `${formatNumericBrushValue(lo)} - ${formatNumericBrushValue(hi)}`;
}

/** Colours for a `useBrushX` brush, to suit the theme. */
export function useBrushStyle(): BrushXStyle {
    const dark = useTheme() === "dark";
    return {
        selectionFill: dark ? "rgba(255, 255, 255, 0.08)" : "rgba(15, 23, 42, 0.08)",
        selectionStroke: dark ? "rgba(229, 231, 235, 0.95)" : "rgba(15, 23, 42, 0.7)",
        handleFill: dark ? "rgba(229, 231, 235, 0.92)" : "rgba(15, 23, 42, 0.88)",
    };
}

/**
 * The brushed range, as a small label over a chart: at the top left by default. Its parent should
 * be `position: relative`; `left` (or `right`) and `top` place it within that.
 */
export default function BrushRangeChip({
    range,
    isDate = false,
    left,
    right,
    top = 4,
}: {
    range: [number, number];
    isDate?: boolean;
    left?: number;
    right?: number;
    top?: number;
}) {
    const dark = useTheme() === "dark";
    return (
        <div
            className="pointer-events-none absolute z-10 rounded px-1 py-0 text-[9px] shadow-sm"
            style={{
                ...(right === undefined ? { left: left ?? 4 } : { right }),
                top,
                border: `1px solid ${dark ? "rgba(71, 85, 105, 0.9)" : "rgba(203, 213, 225, 0.95)"}`,
                background: dark ? "rgba(15, 23, 42, 0.9)" : "rgba(255, 255, 255, 0.92)",
                color: dark ? "rgba(241, 245, 249, 0.95)" : "rgba(15, 23, 42, 0.92)",
            }}
        >
            {formatBrushRange(range[0], range[1], isDate)}
        </div>
    );
}
