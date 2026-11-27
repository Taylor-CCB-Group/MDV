import type { DataType, FieldName } from "@/charts/charts";
import type { FieldSpec, FieldSpecs } from "@/lib/columnTypeHelpers";
import type { ChartArrayCell } from "../hooks/useChartArrayGrid";

/**
 * Chart types that can colour by several columns (`color_fields`) and show a grid of them: any 2D
 * deck scatter plot, and the SpatialData chart.
 */
const COLOR_GRID_CHART_TYPES = new Set([
    "DeckScatter",
    "DeckContourScatter",
    "DeckDensity",
    "wgl_scatter_plot_dev",
    "SpatialDataMdvRegionReact",
]);

export function supportsColorGridMode(chartType: string | undefined) {
    return typeof chartType === "string" && COLOR_GRID_CHART_TYPES.has(chartType);
}

/** How a colour grid makes its cells: one per colour column, or one per category of the first. */
export type ColorGridSplit = "columns" | "categories";

type ColorGridModeConfig = {
    grid?: { by: string; split?: string; layer?: string } | null;
    type?: string;
    dimension?: "2d" | "3d";
};

/**
 * How the chart splits its colour into cells, or null when it draws no colour grid. A grid by
 * `color_fields` with no `split` is a cell per column.
 */
export function colorGridSplit(config: ColorGridModeConfig): ColorGridSplit | null {
    // a chart with no dimension is 2D
    if (!supportsColorGridMode(config.type) || config.dimension === "3d") return null;
    if (config.grid?.by !== "color_fields") return null;
    return config.grid.split === "categories" ? "categories" : "columns";
}

/**
 * The options of a grid that outlast a change to another grid: leaving out empty cells
 * (`hide_empty`), which a grid per region or per category both do.
 */
export function keptGridOptions(previous: object | null | undefined): { hide_empty?: true } {
    return hidesEmptyCells(previous) ? { hide_empty: true } : {};
}

/** Whether a grid leaves out its cells with no rows left after filtering (`hide_empty`). */
export function hidesEmptyCells(grid: object | null | undefined) {
    return !!grid && "hide_empty" in grid && grid.hide_empty === true;
}

/** The `grid` config for a colour grid split this way, keeping the options of the grid it replaces. */
export function colorGridConfig(
    split: ColorGridSplit,
    previous?: object | null,
): { by: "color_fields"; split: ColorGridSplit; hide_empty?: true } {
    return { by: "color_fields", split, ...keptGridOptions(previous) };
}

/**
 * How a grid of a colour column splits it: a cell per category for a column with one category
 * per row, a cell per column (it and those after it) for a numeric one. Null for no column, or a
 * column of another type (`multitext`, `unique`).
 */
export function colorGridSplitFor(column: { datatype: DataType } | undefined): ColorGridSplit | null {
    switch (column?.datatype) {
        case "text":
        case "text16":
            return "categories";
        case "integer":
        case "double":
        case "int32":
            return "columns";
        default:
            return null;
    }
}

/** A grid a chart can show, as its grid settings dropdown and the toolbar offer it. */
export type GridChoice = { name: string; value: string };

/** The grid settings dropdown's choice of no grid. */
export const ONE_VIEW: GridChoice = { name: "One view", value: "none" };

/** The choices a grid settings dropdown can offer for a colour grid. */
export const COLOR_GRID_CHOICES: readonly GridChoice[] = [
    { name: "A cell per colour column", value: "color_fields:columns" },
    { name: "A cell per colour category", value: "color_fields:categories" },
];

/**
 * The colour grid choices to offer: the split the colour column allows (`colorGridSplitFor`), and
 * the chart's current one, so the dropdown can show it.
 */
export function colorGridChoices(available: ColorGridSplit | null, current: ColorGridSplit | null) {
    return COLOR_GRID_CHOICES.filter((choice) => {
        const split = colorGridChoice(choice.value);
        return split === available || split === current;
    });
}

/** The colour grid a grid settings dropdown's value chooses, or null for another kind of grid. */
export function colorGridChoice(value: string): ColorGridSplit | null {
    if (value === "color_fields:columns") return "columns";
    if (value === "color_fields:categories") return "categories";
    return null;
}

/**
 * A grid settings dropdown's value for a chart's grid: what makes the cells, and how for colour
 * (the chart's, or a layer's; see `layerFillGridChoice`).
 */
export function gridChoiceValue(config: ColorGridModeConfig): string {
    const { grid } = config;
    if (!grid) return "none";
    const split = grid.split === "categories" ? "categories" : "columns";
    if (grid.by === "color_fields") return `color_fields:${split}`;
    if (grid.by === "layer_fill" && grid.layer) return `layer_fill:${split}:${grid.layer}`;
    return grid.by;
}

/** The field a legacy `color_by` names: the field itself, or the field of its `column`. */
export function colorByField(colorBy: unknown): FieldName | undefined {
    if (typeof colorBy === "string") return colorBy || undefined;
    if (!colorBy || typeof colorBy !== "object" || !("column" in colorBy)) return undefined;
    const { column } = colorBy;
    if (!column || typeof column !== "object" || !("field" in column)) return undefined;
    return typeof column.field === "string" ? column.field : undefined;
}

/** Some of the rows a layer draws, and where each is in the array they were taken from. */
export type CategoryRows = { rows: Uint32Array; positions: Uint32Array };

/**
 * Splits rows by category, keeping their order. `codes` is a categorical column's data (an index
 * into `values` per row). Categories with no rows are left out.
 */
export function partitionRowsByCategory(
    rows: ArrayLike<number>,
    codes: ArrayLike<number>,
    values: readonly string[],
): Map<string, CategoryRows> {
    const counts = new Uint32Array(values.length);
    for (let i = 0; i < rows.length; i++) {
        const code = codes[rows[i]];
        if (code < values.length) counts[code]++;
    }
    const byCode = Array.from(counts, (count) =>
        count > 0 ? { rows: new Uint32Array(count), positions: new Uint32Array(count) } : null,
    );
    const filled = new Uint32Array(values.length);
    for (let i = 0; i < rows.length; i++) {
        const row = rows[i];
        const code = codes[row];
        const target = byCode[code];
        if (!target) continue;
        const at = filled[code]++;
        target.rows[at] = row;
        target.positions[at] = i;
    }
    const result = new Map<string, CategoryRows>();
    byCode.forEach((target, code) => {
        if (target) result.set(values[code], target);
    });
    return result;
}

/** A colour grid's cell (a category or a column), and the colour of its label. */
export type ColorGridCell = ChartArrayCell & { color: readonly number[] };

const LABEL_FALLBACK_COLOR = [226, 232, 240] as const;

/** A cell per category of a column, in the column's order, labelled in its colour. */
export function colorGridCells(values: readonly string[], colors: readonly unknown[]): ColorGridCell[] {
    return values.map((value, index) => {
        const color = colors[index];
        return {
            key: value,
            label: value,
            color: Array.isArray(color) && color.length >= 3 ? color : LABEL_FALLBACK_COLOR,
        };
    });
}

/**
 * The colour columns with `first` as the colour column. The grid's other columns stay while it is
 * a grid of a cell per column (`keepOthers`); otherwise the colour column is the only one.
 */
export function withColorColumn(specs: FieldSpecs, first: FieldSpec, keepOthers: boolean): FieldSpecs {
    if (!keepOthers) return [first];
    return [first, ...specs.slice(1).filter((spec) => spec !== first)];
}

/** The colour columns with `others` after the colour column, which stays first. */
export function withOtherGridColumns(specs: FieldSpecs, others: FieldSpecs): FieldSpecs {
    const first = specs[0];
    if (first === undefined) return [...others];
    return [first, ...others.filter((spec) => spec !== first)];
}
