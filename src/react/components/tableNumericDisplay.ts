import type { Column, Formatter, MenuCommandItem } from "slickgrid-react";
import { z } from "zod/v4";
import { type TableNumericDisplay, TableNumericDisplaySchema } from "@/charts/schemas/ChartConfigSchema";
import { isDateColumn } from "@/lib/dateFormat";
import { escapeHtml } from "@/utilities/Utilities";
import type { DataColumn, DataType } from "../../charts/charts";

export type ColumnDisplayConfig = Record<string, TableNumericDisplay>;

type ColumnMetadata = Pick<DataColumn<DataType>, "datatype" | "is_date" | "date_unit">;
type ColumnRange = readonly [number, number];

/** Dates are stored as numbers, but neither a bar nor an absolute value means anything for them. */
export function isNumericDisplayColumn(column: ColumnMetadata | null | undefined): boolean {
    if (!column) return false;
    const { datatype } = column;
    return (datatype === "double" || datatype === "integer" || datatype === "int32") && !isDateColumn(column);
}

/**
 * Validate a saved `column_display`. Entries are checked one at a time, so a bad one
 * falls back to a plain number without losing the others.
 */
export function parseColumnDisplay(value: unknown): ColumnDisplayConfig {
    if (value === undefined || value === null) return {};
    const record = z.record(z.string(), z.unknown()).safeParse(value);
    if (!record.success) {
        console.warn("Ignoring invalid table column_display:", record.error.message);
        return {};
    }
    const result: ColumnDisplayConfig = {};
    for (const [field, entry] of Object.entries(record.data)) {
        const parsed = TableNumericDisplaySchema.safeParse(entry);
        if (parsed.success) {
            result[field] = parsed.data;
        } else {
            console.warn(`Ignoring invalid display options for table column '${field}':`, parsed.error.message);
        }
    }
    return result;
}

/** `current` with the display of one column changed. Plain numbers are the default, so are not listed. */
export function withColumnDisplay(
    current: ColumnDisplayConfig | undefined,
    field: string,
    display: TableNumericDisplay,
): ColumnDisplayConfig {
    const next = { ...current };
    if (display.mode === "number") delete next[field];
    else next[field] = display;
    return next;
}

/** The display options offered in a column's header menu; other fixed ranges are set in the chart settings. */
const MENU_CHOICES: { command: string; title: string; display: TableNumericDisplay }[] = [
    { command: "display-number", title: "Show as Number", display: { mode: "number" } },
    { command: "display-bar", title: "Show with Bar", display: { mode: "bar" } },
    { command: "display-diverging", title: "Show with Diverging Bar", display: { mode: "diverging" } },
    {
        command: "display-diverging-unit",
        title: "Show with Diverging Bar (−1 to 1)",
        display: { mode: "diverging", limit: 1 },
    },
];

export function getDisplayForMenuCommand(command: string): TableNumericDisplay | undefined {
    return MENU_CHOICES.find((choice) => choice.command === command)?.display;
}

const getLimit = (display: TableNumericDisplay) => (display.mode === "diverging" ? display.limit : undefined);

/** Header menu items for a numeric column, with a tick on its current display. */
export function getDisplayMenuItems(current: TableNumericDisplay = { mode: "number" }): MenuCommandItem[] {
    const isCurrent = (display: TableNumericDisplay) =>
        display.mode === current.mode && getLimit(display) === getLimit(current);
    // a class without a glyph leaves the icon's space empty: `mdi` alone draws a solid square
    const item = (command: string, title: string, ticked: boolean) => ({
        command,
        title,
        iconCssClass: ticked ? "mdi mdi-check" : "mdv-menu-unticked",
    });
    const items = MENU_CHOICES.map(({ command, title, display }) => item(command, title, isCurrent(display)));
    // a range set in the chart settings is listed too, so that the tick has somewhere to go
    if (!MENU_CHOICES.some(({ display }) => isCurrent(display))) {
        const limit = getLimit(current);
        items.push(item("display-diverging-custom", `Show with Diverging Bar (−${limit} to ${limit})`, true));
    }
    return items;
}

/** The bar of one cell, as percentages of the bar track. */
export type BarSegment = { from: number; to: number; negative: boolean };

const clamp = (value: number, min: number, max: number) => Math.min(Math.max(value, min), max);

export function getBarSegment(
    value: number,
    display: Exclude<TableNumericDisplay, { mode: "number" }>,
    [min, max]: ColumnRange,
): BarSegment {
    if (display.mode === "bar") {
        // a column with a single value has no range: show it as full
        const fraction = max > min ? clamp((value - min) / (max - min), 0, 1) : 1;
        return { from: 0, to: fraction * 100, negative: false };
    }
    const limit = display.limit ?? Math.max(Math.abs(min), Math.abs(max));
    const fraction = limit > 0 ? clamp(value / limit, -1, 1) : 0;
    const end = 50 + fraction * 50;
    return { from: Math.min(50, end), to: Math.max(50, end), negative: value < 0 };
}

/**
 * A short form of a number, so it fits beside its bar: three significant figures, keeping all
 * the integer digits of larger values. The cell's tooltip has the value in full.
 */
export function formatBarValue(value: number): string {
    const magnitude = Math.abs(value);
    if (magnitude >= 1e7) return value.toExponential(2);
    if (Number.isInteger(value)) return String(value);
    if (magnitude < 0.001) return value.toExponential(1);
    const integerDigits = Math.floor(Math.log10(magnitude)) + 1;
    return String(Number(value.toPrecision(Math.max(3, integerDigits))));
}

/**
 * Formatter drawing a number and its bar as a single element: the bar is a CSS gradient
 * (see `.mdv-num-bar` in charts.css), positioned by two custom properties.
 */
export function createBarFormatter(
    display: Exclude<TableNumericDisplay, { mode: "number" }>,
    getRange: () => ColumnRange | undefined,
): Formatter {
    return (_row, _cell, value) => {
        const range = getRange();
        // "missing" and the like have no bar
        if (typeof value !== "number" || !Number.isFinite(value) || !range) {
            return escapeHtml(String(value ?? ""));
        }
        const { from, to, negative } = getBarSegment(value, display, range);
        const diverging = display.mode === "diverging";
        const element = document.createElement("span");
        element.className = `mdv-num-bar${diverging ? " mdv-num-bar-diverging" : ""}${
            negative && diverging ? " mdv-num-bar-negative" : ""
        }`;
        element.style.setProperty("--mdv-bar-from", `${from}%`);
        element.style.setProperty("--mdv-bar-to", `${to}%`);
        element.title = String(value);
        element.textContent = `${diverging && value > 0 ? "+" : ""}${formatBarValue(value)}`;
        return element;
    };
}

/**
 * Set or remove the bar formatter of each grid column to match `display`. Columns are changed
 * in place, so the grid keeps its widths and order.
 *
 * @returns whether any column changed, so the grid needs to redraw its cells.
 */
export function applyColumnDisplay(
    columns: Column[],
    display: ColumnDisplayConfig | undefined,
    getColumn: (field: string) => (ColumnMetadata & { minMax?: ColumnRange }) | undefined,
): boolean {
    let changed = false;
    for (const column of columns) {
        const field = String(column.field);
        const options = display?.[field];
        if (options && options.mode !== "number" && isNumericDisplayColumn(getColumn(field))) {
            column.formatter = createBarFormatter(options, () => getColumn(field)?.minMax);
            changed = true;
        } else if (column.formatter) {
            column.formatter = undefined;
            changed = true;
        }
    }
    return changed;
}
