import type { GuiSpecs } from "@/charts/charts";
import type { ChartGridConfig } from "@/charts/schemas/ChartConfigSchema";
import { g } from "@/lib/utils";
import { densityGridConfig, isGridCrosshairShown, withGridCrosshair } from "./densityGridUtils";

/**
 * The settings folder for a chart that can show its density fields as a grid. The grid's own
 * options are listed only while it is on, as the dialog is built when it opens.
 */
export function getDensityGridSettingsFolder(config: { grid?: ChartGridConfig | null }) {
    const settings: GuiSpecs = [
        g({
            type: "check",
            label: "Show density fields as a grid",
            current_value: config.grid?.by === "densityFields",
            func: (x) => {
                config.grid = densityGridConfig(x, config.grid);
            },
        }),
    ];
    if (config.grid) {
        settings.push(
            g({
                type: "check",
                label: "Mark the pointer in every cell",
                current_value: isGridCrosshairShown(config.grid),
                func: (x) => {
                    config.grid = withGridCrosshair(config.grid, x);
                },
            }),
        );
    }
    return g({ type: "folder", label: "Grid", current_value: settings });
}
