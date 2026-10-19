import type { BaseConfig } from "@/charts/BaseChart";
import type { FieldName } from "@/charts/charts";
import type { AxisConfig2D } from "@/react/scatter_state";
import type { RowOrder } from "./dotPlotLayout";
import type { DotPlotPalette } from "./dotPlotPalettes";

export const DOT_PLOT_REACT_TYPE = "split_dot_plot";

export type DotPlotColorScale = {
    palette: DotPlotPalette;
    log: boolean;
    /**
     * `auto`: from the data, fixed while filtering and recalculated when the rows, fields or
     * split change (as the legacy dot plot does). `manual`: from `min` and `max`.
     */
    mode: "auto" | "manual";
    min?: number;
    max?: number;
};

/**
 * `param` is `[rowColumn, ...fields]`, as for the legacy `dot_plot`: the categories on the y axis,
 * then the fields on the x axis.
 */
export type DotPlotReactConfig = BaseConfig & {
    y_axis_order: RowOrder;
    /** A value counts as expressed when it is greater than this. No setting yet, as in the legacy dot plot. */
    threshold_value?: number;
    /** Category column to split by, with one panel per value in `split_values`; null for one panel. */
    split_by: FieldName | null;
    /** Values of `split_by` to show, in panel order. */
    split_values: string[];
    color_scale: DotPlotColorScale;
    /** `pos` is the legend's saved [left, top] in the chart, as for other charts' color legends. */
    color_legend: { display: boolean; pos?: [number, number] };
    /** `position` is the legend's saved [left, top] in the chart, as in the legacy dot plot. */
    fraction_legend: { display: boolean; position?: [number, number] };
    axis: AxisConfig2D;
};

/** "None" in column dropdowns built with `DataStore.getColumnList(..., true)`. */
export const NO_COLUMN = "__none__";

/**
 * `split_by` and `split_values` from the add-chart dialog's choices: no split for "None", and
 * only values that belong to the chosen column, in the order they were picked.
 */
export function getSplitFromDialog(
    field: unknown,
    values: unknown,
    columnValues: readonly string[] | undefined,
): Pick<DotPlotReactConfig, "split_by" | "split_values"> {
    if (typeof field !== "string" || field === "" || field === NO_COLUMN) {
        return { split_by: null, split_values: [] };
    }
    const known = new Set(columnValues);
    const picked = Array.isArray(values) ? values : [];
    return {
        split_by: field,
        split_values: picked.filter((v): v is string => typeof v === "string" && known.has(v)),
    };
}

/** Fills in defaults, matching the legacy dot plot where it has an equivalent. */
export function adaptDotPlotConfig(original: BaseConfig & Partial<DotPlotReactConfig>): DotPlotReactConfig {
    return {
        ...original,
        y_axis_order: original.y_axis_order ?? "data",
        split_by: original.split_by ?? null,
        split_values: original.split_values ?? [],
        color_scale: {
            ...original.color_scale,
            palette: original.color_scale?.palette ?? "reds",
            log: original.color_scale?.log ?? false,
            mode: original.color_scale?.mode ?? "auto",
        },
        color_legend: { ...original.color_legend, display: original.color_legend?.display ?? true },
        fraction_legend: { ...original.fraction_legend, display: original.fraction_legend?.display ?? true },
        axis: {
            x: { size: 60, tickfont: 10, rotate_labels: true, ...original.axis?.x },
            y: { size: 110, tickfont: 10, rotate_labels: false, ...original.axis?.y },
        },
    };
}

/** Colour range to draw with: the manual one if it's set and valid, otherwise the automatic one. */
export function getColorRange(scale: DotPlotColorScale, autoRange: [number, number] | null): [number, number] | null {
    const { mode, min, max } = scale;
    if (mode === "manual" && typeof min === "number" && typeof max === "number" && Number.isFinite(min + max)) {
        return min <= max ? [min, max] : [max, min];
    }
    return autoRange;
}
