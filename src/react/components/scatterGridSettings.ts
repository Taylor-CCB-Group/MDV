import { g } from "@/lib/utils";
import {
    COLOR_GRID_CHOICES,
    colorGridChoice,
    colorGridConfig,
    gridChoiceValue,
    supportsColorGridMode,
} from "./colorGridUtils";
import { supportsDensityGridMode } from "./densityGridUtils";

/**
 * The deck scatter plot's grid settings: one view, a cell per density field (charts that draw
 * density), a cell per colour column, or a cell per category of the first colour column. Options
 * for the grid's layout belong here as they come.
 */
export function getScatterGridSettingsFolder(config: { type?: string; grid?: { by: string; split?: string } | null }) {
    const choices = [
        { name: "One view", value: "none" },
        ...(supportsDensityGridMode(config.type) ? [{ name: "A cell per density field", value: "densityFields" }] : []),
        ...(supportsColorGridMode(config.type) ? COLOR_GRID_CHOICES : []),
    ];
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
                values: [choices, "name", "value"],
                func: (value: string) => {
                    const colorSplit = colorGridChoice(value);
                    if (colorSplit) config.grid = colorGridConfig(colorSplit);
                    else config.grid = value === "densityFields" ? { by: value } : null;
                },
            }),
        ],
    });
}
