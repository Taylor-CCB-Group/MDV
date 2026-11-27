import { autorun, observable, runInAction } from "mobx";
import type { Column } from "slickgrid-react";
import { describe, expect, test, vi } from "vitest";

vi.mock("@/react/components/TableChartReactComponent", () => ({
    default: () => null,
}));
vi.mock("@/react/components/BaseReactChart", () => ({
    BaseReactChart: class {},
}));
vi.mock("@/charts/BaseChart", () => ({
    // biome-ignore lint/complexity/noStaticOnlyClass: stands in for BaseChart's static registry
    default: class {
        static types: Record<string, unknown> = {};
    },
}));

import { getChartConfigSchema } from "@/charts/schemas/ChartConfigRegistry";
import { TableChartReactConfigSchema } from "@/charts/schemas/ChartConfigSchema";
import { adaptConfig, type TableChartReactConfig } from "@/react/components/TableChartReactWrapper";
import {
    applyColumnDisplay,
    formatBarValue,
    getBarSegment,
    parseColumnDisplay,
} from "@/react/components/tableNumericDisplay";
import { getNumericDisplaySettings } from "@/react/components/tableNumericDisplaySettings";

const savedConfig = (extra: Record<string, unknown> = {}): TableChartReactConfig => ({
    id: "table",
    type: "table_chart_react",
    title: "Features",
    legend: "",
    size: [400, 300],
    param: ["name", "rho"],
    ...extra,
});

describe("table config", () => {
    test("is registered", () => {
        expect(getChartConfigSchema("table_chart_react")).toBe(TableChartReactConfigSchema);
    });

    test("accepts a view saved without display options", () => {
        const config = savedConfig({
            include_index: false,
            column_widths: { rho: 150 },
            order: { name: 0, rho: 1 },
            sort: { columnId: "rho", ascending: false },
        });

        expect(TableChartReactConfigSchema.safeParse(config).success).toBe(true);
        expect(adaptConfig(config)).toMatchObject({
            sort: { columnId: "rho", ascending: false },
            column_display: {},
        });
    });

    test("accepts display options and an absolute sort", () => {
        const column_display = {
            rho: { mode: "diverging", limit: 1 },
            area: { mode: "bar" },
            count: { mode: "number" },
        };
        const sort = { columnId: "rho", ascending: false, absolute: true };
        const config = savedConfig({ column_display, sort });

        expect(TableChartReactConfigSchema.safeParse(config).success).toBe(true);
        expect(adaptConfig(config)).toMatchObject({ column_display, sort });
    });

    test.each([
        { mode: "sparkline" },
        { mode: "diverging", limit: 0 },
        { mode: "diverging", limit: "1" },
        "bar",
        null,
    ])("rejects the display option %j, keeping the other columns", (invalid) => {
        vi.spyOn(console, "warn").mockImplementation(() => {});
        const column_display = { bad: invalid, rho: { mode: "bar" } };

        expect(TableChartReactConfigSchema.safeParse(savedConfig({ column_display })).success).toBe(false);
        expect(parseColumnDisplay(column_display)).toEqual({ rho: { mode: "bar" } });
        vi.restoreAllMocks();
    });

    test.each([undefined, null, "bar", [{ mode: "bar" }]])("parses a column_display of %j as empty", (value) => {
        vi.spyOn(console, "warn").mockImplementation(() => {});
        expect(parseColumnDisplay(value)).toEqual({});
        vi.restoreAllMocks();
    });

    test("drops an absolute flag that is not a boolean", () => {
        const config = savedConfig({ sort: { columnId: "rho", ascending: true, absolute: "yes" } });

        expect(adaptConfig(config).sort).toEqual({ columnId: "rho", ascending: true });
    });
});

describe("getBarSegment", () => {
    test("scales a bar to the column's range", () => {
        const display = { mode: "bar" } as const;
        expect(getBarSegment(10, display, [10, 30])).toEqual({ from: 0, to: 0, negative: false });
        expect(getBarSegment(25, display, [10, 30])).toEqual({ from: 0, to: 75, negative: false });
        // values edited after the range was measured stay inside the track
        expect(getBarSegment(40, display, [10, 30]).to).toBe(100);
        expect(getBarSegment(5, display, [5, 5]).to).toBe(100);
    });

    test("centres a diverging bar on zero, scaled to the larger side of the data", () => {
        const display = { mode: "diverging" } as const;
        expect(getBarSegment(0.4, display, [-0.8, 0.4])).toEqual({ from: 50, to: 75, negative: false });
        expect(getBarSegment(-0.8, display, [-0.8, 0.4])).toEqual({ from: 0, to: 50, negative: true });
        expect(getBarSegment(0, display, [0, 0])).toEqual({ from: 50, to: 50, negative: false });
    });

    test("uses a fixed symmetric range when given", () => {
        const display = { mode: "diverging", limit: 1 } as const;
        expect(getBarSegment(-0.5, display, [-0.6, 0.6])).toEqual({ from: 25, to: 50, negative: true });
        expect(getBarSegment(3, display, [-0.6, 3])).toEqual({ from: 50, to: 100, negative: false });
    });
});

describe("formatBarValue", () => {
    test.each([
        [0.46000000834465027, "0.46"],
        [-0.123456, "-0.123"],
        [0.001234, "0.00123"],
        [0.00001234, "1.2e-5"],
        [12345.678, "12346"],
        [123456, "123456"],
        [123456789, "1.23e+8"],
        [0, "0"],
    ])("shortens %d to %s", (value, expected) => {
        expect(formatBarValue(value)).toBe(expected);
    });
});

describe("applyColumnDisplay", () => {
    const columnIndex: Record<string, any> = {
        name: { datatype: "text" },
        rho: { datatype: "double", minMax: [-1, 1] },
        visit: { datatype: "int32", is_date: true, minMax: [0, 100] },
    };
    const getColumn = (field: string) => columnIndex[field];
    const makeColumns = (): Column[] => ["name", "rho", "visit"].map((field) => ({ id: field, field }));

    const render = (column: Column | undefined, value: unknown) =>
        column?.formatter?.(0, 0, value, column, {}, {} as never);

    test("leaves columns alone when there are no display options", () => {
        const columns = makeColumns();
        expect(applyColumnDisplay(columns, {}, getColumn)).toBe(false);
        expect(columns).toEqual(makeColumns());
    });

    test("draws bars only for numeric columns, and removes them again", () => {
        const columns = makeColumns();
        const bar = { mode: "bar" } as const;

        expect(applyColumnDisplay(columns, { name: bar, rho: bar, visit: bar }, getColumn)).toBe(true);
        expect(columns.map((column) => Boolean(column.formatter))).toEqual([false, true, false]);

        expect(applyColumnDisplay(columns, { rho: { mode: "number" } }, getColumn)).toBe(true);
        expect(columns).toEqual(makeColumns());
    });

    test("renders a signed value and its bar as one element", () => {
        const columns = makeColumns();
        applyColumnDisplay(columns, { rho: { mode: "diverging", limit: 1 } }, getColumn);

        const positive = render(columns[1], 0.5) as HTMLElement;
        expect(positive.textContent).toBe("+0.5");
        expect(positive.childElementCount).toBe(0);
        expect(positive.style.getPropertyValue("--mdv-bar-from")).toBe("50%");
        expect(positive.style.getPropertyValue("--mdv-bar-to")).toBe("75%");
        expect(positive.classList.contains("mdv-num-bar-negative")).toBe(false);

        const negative = render(columns[1], -0.5) as HTMLElement;
        expect(negative.textContent).toBe("-0.5");
        expect(negative.classList.contains("mdv-num-bar-negative")).toBe(true);

        expect(render(columns[1], "missing")).toBe("missing");
    });
});

describe("numeric display settings", () => {
    const columnIndex: Record<string, any> = {
        name: { name: "Name", datatype: "text" },
        rho: { name: "Spearman rho", datatype: "double" },
        n: { name: "n", datatype: "integer" },
    };

    function controls(folder: ReturnType<typeof getNumericDisplaySettings>, index = 0) {
        const column = folder.current_value[index];
        if (column?.type !== "folder") throw new Error("expected a folder per column");
        const [mode, range] = column.current_value;
        if (mode?.type !== "radiobuttons" || range?.type !== "text") throw new Error("unexpected controls");
        return { mode, range };
    }

    /** As the dialog's text box does: set the value, then call `func`, in one action. */
    function type(range: ReturnType<typeof controls>["range"], text: string) {
        runInAction(() => {
            range.current_value = text;
            range.func?.(text);
        });
    }

    test("has a folder per numeric column, following the table's columns", () => {
        const config = observable({ param: ["name", "rho"] });
        const folder = getNumericDisplaySettings(config, columnIndex);
        let labels: string[] = [];
        const dispose = autorun(() => {
            labels = folder.current_value.map((child) => child.label);
        });
        expect(labels).toEqual(["Spearman rho"]);

        runInAction(() => {
            config.param = ["name", "rho", "n"];
        });
        expect(labels).toEqual(["Spearman rho", "n"]);
        dispose();
    });

    test("writes the chosen mode and range to the config", () => {
        const config = observable({ param: ["rho"], column_display: {} });
        const folder = getNumericDisplaySettings(config, columnIndex);
        const { mode, range } = controls(folder);

        mode.func?.("diverging");
        expect(mode.current_value).toBe("diverging");
        expect(config.column_display).toEqual({ rho: { mode: "diverging" } });

        type(range, "1");
        expect(config.column_display).toEqual({ rho: { mode: "diverging", limit: 1 } });

        // partial input is left in the box and not saved
        type(range, "0.");
        expect(range.current_value).toBe("0.");
        expect(config.column_display).toEqual({ rho: { mode: "diverging", limit: 1 } });

        // text that means the saved range is left as typed
        type(range, "2.0");
        expect(range.current_value).toBe("2.0");
        expect(config.column_display).toEqual({ rho: { mode: "diverging", limit: 2 } });

        // blank goes back to the data range
        type(range, "");
        expect(range.current_value).toBe("");
        expect(config.column_display).toEqual({ rho: { mode: "diverging" } });

        mode.func?.("number");
        expect(config.column_display).toEqual({});
    });

    test("follows a display changed elsewhere, such as from the header menu", () => {
        const config = observable({ param: ["rho"], column_display: {} });
        const folder = getNumericDisplaySettings(config, columnIndex);
        const { mode, range } = controls(folder);
        let shown: string[] = [];
        const dispose = autorun(() => {
            shown = [mode.current_value, range.current_value];
        });

        runInAction(() => {
            config.column_display = { rho: { mode: "diverging", limit: 1 } };
        });
        expect(shown).toEqual(["diverging", "1"]);

        // switching away and back keeps that range
        mode.func?.("bar");
        mode.func?.("diverging");
        expect(config.column_display).toEqual({ rho: { mode: "diverging", limit: 1 } });

        runInAction(() => {
            config.column_display = { rho: { mode: "diverging" } };
        });
        expect(shown).toEqual(["diverging", ""]);
        dispose();
    });
});
