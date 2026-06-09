import type DataStore from "@/datastore/DataStore";
import { isColumnNumeric, isColumnText } from "@/utilities/Utilities";
import type {
    ColorLegendBuildConfig,
    ColorLegendCategoricalItem,
    ColorLegendSpec,
} from "./types";

/**
 * Builds a serializable color legend spec from a datastore column.
 * Replaces DataStore.getColorLegend DOM construction for the color_by legend path.
 */
export function buildColorLegendSpec(
    dataStore: DataStore,
    column: string,
    config: ColorLegendBuildConfig = {},
): ColorLegendSpec | null {
    const c = dataStore.columnIndex[column];
    if (!c) {
        return null;
    }

    const overrideValues = config.overrideValues ?? config.overideValues;
    const colorConfig = overrideValues
        ? { ...config, overideValues: overrideValues }
        : config;
    const rawColors = dataStore.getColumnColors(column, colorConfig);
    if (!rawColors) {
        return null;
    }
    const colors = rawColors.map(color => String(color));

    const name = config.name ?? c.name;

    if (isColumnNumeric(c)) {
        const [min, max] = dataStore.getMinMaxForColumn(column);
        let range: [number, number] = [min, max];
        if (overrideValues) {
            const ov = overrideValues;
            range = [
                ov.min == null ? min : ov.min,
                ov.max == null ? max : ov.max,
            ];
        }
        return {
            kind: "continuous",
            label: name,
            colors,
            range,
        };
    }

    if (isColumnText(c)) {
        const values = c.values ?? [];
        const items: ColorLegendCategoricalItem[] = values.map(
            (value, i) => ({
                color: String(colors[i]),
                name: String(value),
                value: String(value),
            }),
        );
        return {
            kind: "categorical",
            label: name,
            column,
            items,
        };
    }

    return null;
}
