import { comparer, reaction, runInAction } from "mobx";
import { useEffect, useMemo } from "react";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig, ColorConfig } from "@/charts/BaseChart";
import type { DataColumn, DataType, FieldName, GridBySpec } from "@/charts/charts";
import { type FieldSpec, type FieldSpecs, flattenFields } from "@/lib/columnTypeHelpers";
import { isDateColumn } from "@/lib/dateFormat";
import { g } from "@/lib/utils";
import {
    type ColorGridSplit,
    colorByField,
    colorGridConfig,
    colorGridSplit,
    colorGridSplitFor,
    supportsColorGridMode,
} from "./components/colorGridUtils";
import { useChart } from "./context";
import { useConfig, useFieldSpecs } from "./hooks";
import { buildColorLegendSpec } from "./legend/color_legend/buildColorLegendSpec";
import type { ColorLegendSpec } from "./legend/color_legend/types";

/**
 * Colour by several columns, for React charts that draw points (the deck scatter plot and the
 * SpatialData chart). It replaces `color_by` there, and BaseChart's handling of it: these charts
 * load the columns themselves, as they do `densityFields`, so a live query (an active link) is
 * just a value of `color_fields`. One view is coloured by the first column; a grid
 * (`grid: { by: "color_fields" }`) can draw a cell per column, or a cell per category of the first.
 */
export type ColorFieldsConfig = Omit<ColorConfig, "color_by"> & {
    color_fields: FieldSpecs;
    grid?: { by: string; split?: string } | null;
    type?: string;
    dimension?: "2d" | "3d";
};

/** A row's colour, as an rgb array. */
export type RowColor = (row: number) => [r: number, g: number, b: number];

/**
 * Gives a config saved with `color_by` (any chart saved before `color_fields` existed, or converted
 * from another type) its `color_fields`, and leaves `color_by` out: these charts no longer read or
 * write it. A config with `color_fields` keeps them.
 */
export function migrateColorBy<T extends { color_by?: unknown; color_fields?: FieldSpecs }>(
    config: T,
): T & { color_fields: FieldSpecs } {
    const field = colorByField(config.color_by);
    const color_fields = config.color_fields ?? (field ? [field] : []);
    delete config.color_by;
    return Object.assign(config, { color_fields });
}

/** The columns `color_fields` names now, active links resolved; an empty list when there are none. */
export function colorFieldNames(colorFields: FieldSpecs | undefined): FieldName[] {
    return colorFields && colorFields.length > 0 ? flattenFields(colorFields) : [];
}

/**
 * The column the colour legend shows: the first colour column, unless each cell of a grid has a
 * column of its own (each cell is labelled with it instead).
 */
export function colorLegendField(config: ColorFieldsConfig): FieldName | undefined {
    if (colorGridSplit(config) === "columns") return undefined;
    return colorFieldNames(config.color_fields)[0];
}

/**
 * The colour function for a column, with the chart's colour settings (log scale, trim, zero as
 * missing, hide missing). The column must be loaded.
 */
export function colorFunctionFor(chart: BaseChart<BaseConfig>, field: FieldName): RowColor {
    const conf = {
        asArray: true,
        overideValues: {
            colorLogScale: chart.config.log_color_scale,
            fallbackOnZero: chart.config.fallbackOnZero,
            hideMissing: chart.config.hideMissing,
        },
    };
    chart._addTrimmedColor(field, conf);
    return chart.dataStore.getColorFunction(field, conf);
}

export type ColorFields = {
    /** The colour columns, loaded, in the order of `color_fields`. */
    columns: DataColumn<DataType>[];
    /** Each column's colour function, aligned with `columns`. */
    colorFunctions: RowColor[];
};

const NO_FIELDS: FieldSpecs = [];

/**
 * The chart's colour columns and their colour functions. A chart without `color_fields` (one that
 * still uses `color_by`) has none. When a colour column's colours change in place, the chart
 * replaces `color_fields` with a copy, which loads the columns again and so remakes these.
 */
export function useColorFields(): ColorFields {
    const chart = useChart();
    const config = useConfig<Partial<ColorFieldsConfig>>();
    const columns = useFieldSpecs(config.color_fields ?? NO_FIELDS);
    const { log_color_scale, trim_color_scale, fallbackOnZero, hideMissing } = config;
    // biome-ignore lint/correctness/useExhaustiveDependencies: the settings are read through `chart`
    const colorFunctions = useMemo(
        () => columns.map((column) => colorFunctionFor(chart, column.field)),
        [chart, columns, log_color_scale, trim_color_scale, fallbackOnZero, hideMissing],
    );
    return { columns, colorFunctions };
}

/**
 * Redraws the colour legend once its column has loaded. The chart redraws it when the legend
 * column changes (see `setUpColorFields`), but a column an active link has just named is loaded
 * after that, and the legend is built from the loaded column. Call once per chart, from its root.
 */
export function useColorFieldsLegend() {
    const chart = useChart();
    const config = useConfig<ColorFieldsConfig>();
    const { columns } = useColorFields();
    const field = colorLegendField(config);
    const column = columns.find((c) => c.field === field);
    useEffect(() => {
        if (column) chart.setColorLegend();
    }, [chart, column]);
}

/**
 * The colour of the chart's points in one view, and the column it comes from: the first colour
 * column for a chart with `color_fields`, else the chart's `colorBy` (set from `color_by`).
 */
export function useChartColorBy(): { colorBy: RowColor | undefined; colorColumn: DataColumn<DataType> | undefined } {
    // the charts that still use `color_by` set `colorBy` from it
    const chart = useChart<BaseConfig, BaseChart<BaseConfig> & { colorBy?: RowColor | null }>();
    const config = useConfig<Partial<ColorFieldsConfig> & { color_by?: unknown }>();
    const { columns, colorFunctions } = useColorFields();
    const legacyColumns = useFieldSpecs(config.color_fields ? undefined : colorByField(config.color_by));
    if (config.color_fields) return { colorBy: colorFunctions[0], colorColumn: columns[0] };
    return { colorBy: chart.colorBy ?? undefined, colorColumn: legacyColumns[0] };
}

type TrimChoice = NonNullable<ColorConfig["trim_color_scale"]>;
const TRIMS: readonly TrimChoice[] = ["none", "0.001", "0.01", "0.05"];
const isTrim = (value: string): value is TrimChoice => TRIMS.some((trim) => trim === value);
const TRIM_CHOICES: [string, string][] = TRIMS.map((trim) => [trim === "none" ? "No Trim" : trim, trim]);

/** Where a chart looks up its columns' types: its datastore. */
export type ColumnTypes = { columnIndex: Record<string, { datatype: DataType } | undefined> };

/**
 * How a grid of the chart's colour would be split now, by the type of its first colour column
 * (see `colorGridSplitFor`), or null when it cannot be a grid of it.
 */
export function availableColorGridSplit(
    config: { color_fields?: FieldSpecs },
    columns: ColumnTypes,
): ColorGridSplit | null {
    const first = colorFieldNames(config.color_fields)[0];
    return colorGridSplitFor(first ? columns.columnIndex[first] : undefined);
}

/** For a chart that can be a grid of its colour: turns the grid on (split by the column's type) and off. */
function colorFieldsGridBy(config: ColorFieldsConfig, columns: ColumnTypes): GridBySpec | undefined {
    if (!supportsColorGridMode(config.type)) return undefined;
    return {
        get: () => colorGridSplit(config),
        available: () => availableColorGridSplit(config, columns),
        set: (on) => {
            const split = on ? availableColorGridSplit(config, columns) : null;
            if (split) config.grid = colorGridConfig(split);
            else if (config.grid?.by === "color_fields") config.grid = null;
        },
    };
}

/** The colour settings of a chart with `color_fields`, in place of BaseChart's "Color Settings". */
export function getColorFieldsSettingsFolder(config: ColorFieldsConfig, columns: ColumnTypes) {
    return g({
        type: "folder",
        label: "Color Settings",
        current_value: [
            g({
                type: "multicolumn",
                label: "Color By",
                current_value: config.color_fields,
                // cannot colour by unique (at the moment)
                columnType: ["int32", "text", "integer", "double", "text16", "multitext"],
                func: (x: FieldSpecs) => {
                    config.color_fields = x;
                },
                gridBy: colorFieldsGridBy(config, columns),
            }),
            g({
                type: "check",
                label: "Show Color Legend",
                current_value: config.color_legend?.display ?? true,
                func: (x: boolean) => {
                    config.color_legend = { ...config.color_legend, display: x };
                },
            }),
            g({
                type: "check",
                label: "SymLog Color Scale",
                current_value: config.log_color_scale ?? false,
                func: (x: boolean) => {
                    config.log_color_scale = x;
                },
            }),
            g({
                type: "check",
                label: "Treat zero as missing",
                current_value: config.fallbackOnZero ?? false,
                func: (x: boolean) => {
                    config.fallbackOnZero = x;
                },
            }),
            g({
                type: "check",
                label: "Hide missing values",
                current_value: config.hideMissing ?? false,
                func: (x: boolean) => {
                    config.hideMissing = x;
                },
            }),
            g({
                type: "radiobuttons",
                label: "Trim Color Scale to Percentile",
                current_value: config.trim_color_scale ?? "none",
                choices: TRIM_CHOICES,
                func: (v) => {
                    if (isTrim(v)) config.trim_color_scale = v;
                },
            }),
        ],
    });
}

/** The config of a chart that colours its points by `color_fields`. */
type ChartColorFieldsConfig = BaseConfig & ColorFieldsConfig;

/**
 * The colour legend of a chart with `color_fields`: for the column `colorLegendField` names, with
 * the chart's log scale and trim, as BaseChart makes it for `color_by`.
 */
export function colorFieldsLegendSpec<T extends ChartColorFieldsConfig>(chart: BaseChart<T>): ColorLegendSpec | null {
    const config: ChartColorFieldsConfig = chart.config;
    const field = colorLegendField(config);
    const column = field ? chart.dataStore.columnIndex[field] : undefined;
    if (!field || !column) return null;
    const conf = {
        overideValues: {
            // Log color scales are not meaningful for calendar day numbers.
            colorLogScale: isDateColumn(column) ? false : config.log_color_scale,
        },
    };
    chart._addTrimmedColor(field, conf);
    return buildColorLegendSpec(chart.dataStore, field, conf);
}

/**
 * Sets up a chart that colours its points by `color_fields`, once BaseChart's constructor has run:
 *
 * - `colorByColumn`, which links and other callers use to colour a chart by one column, sets
 *   `color_fields` to that column. It is set on the chart rather than declared on its class, so
 *   that BaseChart does not wrap it to load the column and track queries: the chart's components
 *   load `color_fields` themselves, live queries included.
 * - the colour legend follows the legend column and the settings that change it.
 * - a grid of the colour keeps the split its first column's type allows, when the columns change
 *   other than by its picker (an active link, `colorByColumn`).
 */
export function setUpColorFields<T extends ChartColorFieldsConfig>(chart: BaseChart<T>) {
    const config: ChartColorFieldsConfig = chart.config;
    chart.colorByColumn = (column) => {
        runInAction(() => {
            config.color_fields = column ? [column] : [];
            // a column pointer link writes the property it targets before calling this
            config.color_by = undefined;
        });
    };
    chart.reactionDisposers.push(
        reaction(
            () => [
                colorLegendField(config),
                config.log_color_scale,
                config.trim_color_scale,
                config.color_legend?.display,
            ],
            () => chart.setColorLegend(),
            { equals: comparer.structural },
        ),
        reaction(
            () => (colorGridSplit(config) ? availableColorGridSplit(config, chart.dataStore) : null),
            (split) => {
                // none (a link not yet resolved, or a column of another type) leaves the grid be
                if (split && colorGridSplit(config) !== split) config.grid = colorGridConfig(split);
            },
        ),
    );
}

/**
 * After a column's data or colours change: the colour columns are loaded again (by replacing
 * `color_fields` with a copy), which remakes their colour functions, and the legend is redrawn.
 */
export function onColorFieldsDataChanged<T extends ChartColorFieldsConfig>(
    chart: BaseChart<T>,
    columns: readonly FieldName[],
) {
    const config: ChartColorFieldsConfig = chart.config;
    const fields = colorFieldNames(config.color_fields);
    if (!fields.some((field) => columns.includes(field))) return;
    runInAction(() => {
        config.color_fields = config.color_fields.slice();
    });
    chart.setColorLegend();
}

/** Stops colouring by a column that has been removed from the datastore. */
export function onColorFieldRemoved<T extends ChartColorFieldsConfig>(chart: BaseChart<T>, removed: FieldName) {
    const config: ChartColorFieldsConfig = chart.config;
    const { color_fields } = config;
    if (!color_fields.some((spec: FieldSpec) => spec === removed)) return;
    runInAction(() => {
        config.color_fields = color_fields.filter((spec: FieldSpec) => spec !== removed);
    });
}
