/**
 * MDV's colour scheme for one column, expressed as viewer config.
 *
 * The viewer resolves a colour column itself — it reads the SpatialData table, so
 * it needs nothing from us to do it. What it cannot know is which colours MDV has
 * already chosen for that column, and a layer that draws `cell_type` in a
 * different palette from every other chart over the same data is worse than one
 * that does not draw it at all.
 *
 * So MDV supplies the SCHEME and the viewer does the work. Everything here is
 * plain JSON that travels in the saved render stack, which means a stack saved
 * from MDV opens in any viewer with MDV's colours intact — where the old approach
 * (MDV computing every feature's colour into `featureState`) produced a stack
 * that only looked right inside MDV, and cost a `Record<string, [r,g,b,a]>` with
 * one entry per cell to produce it.
 *
 * Columns that are NOT in the table's obs — linked gene scores, `mdv_cell_id`,
 * anything computed at runtime — cannot go down this path at all: the viewer has
 * no way to read them. Those keep the per-feature route. See `table_association`.
 */

import type { LayerConfig } from "@spatialdata/vis";

import type DataStore from "@/datastore/DataStore";
import { isDatatypeNumeric } from "@/lib/utils";

type ShapesLayerConfig = Extract<LayerConfig, { type: "shapes" }>;
export type FillColorByColumn = NonNullable<ShapesLayerConfig["fillColorByColumn"]>;

type RgbColor = [number, number, number];

const CATEGORICAL_DATATYPES = new Set(["text", "text16", "multitext"]);

function isRgb(value: unknown): value is RgbColor {
    return (
        Array.isArray(value) &&
        value.length >= 3 &&
        typeof value[0] === "number" &&
        typeof value[1] === "number" &&
        typeof value[2] === "number"
    );
}

function rgbList(colors: unknown): RgbColor[] {
    if (!Array.isArray(colors)) return [];
    return colors.flatMap((color) => (isRgb(color) ? [[color[0], color[1], color[2]] as RgbColor] : []));
}

/**
 * How MDV colours `columnName`, as a `fillColorByColumn` the viewer can execute.
 *
 * `undefined` means "MDV has no opinion it can express" — an unknown column, an
 * unsupported datatype, or a numeric column whose min/max is not known yet
 * (`minMax` is only computed once the data lands). The caller should fall back to
 * the viewer's own defaults rather than treat that as an error; the colours will
 * simply be the viewer's rather than MDV's until the column loads.
 */
export function fillColorSchemeFromDataStore(dataStore: DataStore, columnName: string): FillColorByColumn | undefined {
    const column = dataStore.columnIndex[columnName];
    if (!column) return undefined;

    if (CATEGORICAL_DATATYPES.has(column.datatype)) {
        const values: string[] = Array.isArray(column.values) ? column.values : [];
        if (values.length === 0) return undefined;
        const colors = rgbList(dataStore.getColumnColors(columnName, { asArray: true }));
        const byValue: Record<string, RgbColor> = {};
        for (const [index, value] of values.entries()) {
            const color = colors[index];
            // Values are matched as strings, because that is the canonical form the
            // viewer normalises a cell to — a category spelled `3` in an integer-coded
            // column is `"3"` on both sides.
            if (color) byValue[String(value)] = color;
        }
        if (Object.keys(byValue).length === 0) return undefined;
        // `mode` states what the palette IS, rather than asking the viewer to
        // re-derive it: MDV knows this column is categorical because it read the
        // datatype, and a `byValue` palette is meaningless under a continuous mode.
        return { columnName, mode: "categorical", categoricalPalette: { byValue } };
    }

    if (isDatatypeNumeric(column.datatype)) {
        const minMax = column.minMax;
        if (!Array.isArray(minMax) || minMax.length < 2) return undefined;
        // `getColumnColors` returns the ramp already sampled into evenly spaced bins
        // across [min, max] — and already remapped if the column is on a log scale.
        // Handing those over as the ramp stops reproduces MDV's scale exactly, which
        // is why nothing here sets `numericScale`: doing both would apply the log
        // twice. The viewer interpolates between the bins rather than stepping, so
        // the result is MDV's ramp, smoothed.
        const stops = rgbList(dataStore.getColumnColors(columnName, { asArray: true }));
        if (stops.length < 2) return undefined;
        return {
            columnName,
            mode: "continuous",
            numericRamp: stops as [RgbColor, RgbColor, ...RgbColor[]],
            numericDomain: [minMax[0], minMax[1]],
        };
    }

    return undefined;
}
