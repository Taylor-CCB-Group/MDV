import type { FieldSpecs } from "@/lib/columnTypeHelpers";
import { g } from "@/lib/utils";
import { availableColorGridSplit, type ColumnTypes } from "../colorFields";
import {
    colorGridChoice,
    colorGridChoices,
    colorGridConfig,
    colorGridSplit,
    type GridChoice,
    gridChoiceValue,
    ONE_VIEW,
    supportsColorGridMode,
} from "./colorGridUtils";
import { supportsDensityGridMode } from "./densityGridUtils";

/**
 * In a grid with a cell per category: whether each cell also shows the other categories, dimmed
 * (points in grey; a layer's shapes or labels faded).
 */
export function getGreyOtherCategoriesSetting(config: { grey_other_categories?: boolean | null }) {
    return g({
        type: "check",
        label: "A cell per category: also show the other categories, dimmed",
        current_value: !!config.grey_other_categories,
        func: (x: boolean) => {
            config.grey_other_categories = x;
        },
    });
}

/**
 * In a grid with a cell per category or per region: whether cells with no rows left after
 * filtering are left out (`grid.hide_empty`), which outlasts a change between such grids.
 */
export function getHideEmptyCellsSetting(config: { grid?: { by: string; hide_empty?: boolean } | null }) {
    return g({
        type: "check",
        label: "A cell per category or region: leave out those with no rows after filtering",
        current_value: config.grid?.hide_empty ?? false,
        func: (x: boolean) => {
            if (config.grid) config.grid = { ...config.grid, hide_empty: x || undefined };
        },
    });
}

type ScatterGridSettingsConfig = {
    type?: string;
    grid?: { by: string; split?: string; hide_empty?: boolean } | null;
    color_fields?: FieldSpecs;
    grey_other_categories?: boolean | null;
};

/**
 * The grids the deck scatter plot can show, besides one view: a cell per density field (charts
 * that draw density), or a grid of the colour, as the colour column's type allows.
 */
export function scatterGridChoices(config: ScatterGridSettingsConfig, columns: ColumnTypes): GridChoice[] {
    return [
        ...(supportsDensityGridMode(config.type) ? [{ name: "A cell per density field", value: "densityFields" }] : []),
        ...(supportsColorGridMode(config.type)
            ? colorGridChoices(availableColorGridSplit(config, columns), colorGridSplit(config))
            : []),
    ];
}

/** Shows the grid a choice names (`scatterGridChoices`), or one view for "none". */
export function setScatterGrid(config: ScatterGridSettingsConfig, value: string) {
    const colorSplit = colorGridChoice(value);
    if (colorSplit) config.grid = colorGridConfig(colorSplit, config.grid);
    else config.grid = value === "densityFields" ? { by: value } : null;
}

/**
 * The deck scatter plot's grid settings: one view, a cell per density field (charts that draw
 * density), or a grid of the colour: a cell per colour column, or per category of a categorical
 * colour column. Options for the grid's layout belong here as they come.
 */
export function getScatterGridSettingsFolder(config: ScatterGridSettingsConfig, columns: ColumnTypes) {
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
                get values(): [GridChoice[], string, string] {
                    return [[ONE_VIEW, ...scatterGridChoices(config, columns)], "name", "value"];
                },
                set values(_values) {
                    // made from the config
                },
                func: (value: string) => setScatterGrid(config, value),
            }),
            getGreyOtherCategoriesSetting(config),
            getHideEmptyCellsSetting(config),
        ],
    });
}
