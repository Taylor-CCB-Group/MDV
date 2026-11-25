import type { FieldSpecs } from "@/lib/columnTypeHelpers";
import { g } from "@/lib/utils";
import { availableColorGridSplit, type ColumnTypes } from "../colorFields";
import {
    colorGridChoice,
    colorGridChoices,
    colorGridConfig,
    colorGridSplit,
    gridChoiceValue,
    supportsColorGridMode,
} from "./colorGridUtils";
import { supportsDensityGridMode } from "./densityGridUtils";

/**
 * The deck scatter plot's grid settings: one view, a cell per density field (charts that draw
 * density), or a grid of the colour: a cell per colour column, or per category of a categorical
 * colour column. Options for the grid's layout belong here as they come.
 */
export function getScatterGridSettingsFolder(
    config: { type?: string; grid?: { by: string; split?: string } | null; color_fields?: FieldSpecs },
    columns: ColumnTypes,
) {
    return g({
        type: "folder",
        label: "Grid",
        current_value: [
            g({
                type: "dropdown",
                label: "Show",
                // read from the config, as the colour pickers' "Grid by this" also sets it
                get current_value() {
                    return gridChoiceValue(config);
                },
                set current_value(_value: string) {
                    // `func` sets the grid, which the getter reads back
                },
                // the colour choice follows the colour column's type
                get values(): [{ name: string; value: string }[], string, string] {
                    const choices = [
                        { name: "One view", value: "none" },
                        ...(supportsDensityGridMode(config.type)
                            ? [{ name: "A cell per density field", value: "densityFields" }]
                            : []),
                        ...(supportsColorGridMode(config.type)
                            ? colorGridChoices(availableColorGridSplit(config, columns), colorGridSplit(config))
                            : []),
                    ];
                    return [choices, "name", "value"];
                },
                set values(_values) {
                    // made from the config
                },
                func: (value: string) => {
                    const colorSplit = colorGridChoice(value);
                    if (colorSplit) config.grid = colorGridConfig(colorSplit);
                    else config.grid = value === "densityFields" ? { by: value } : null;
                },
            }),
        ],
    });
}
