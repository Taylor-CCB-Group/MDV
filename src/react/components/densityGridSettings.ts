import { g } from "@/lib/utils";
import { densityGridConfig } from "./densityGridUtils";

/**
 * The settings folder for a chart that can show its density fields as a grid. Options for the
 * grid's layout belong here as they come.
 */
export function getDensityGridSettingsFolder(config: { grid?: { by: string } | null }) {
    return g({
        type: "folder",
        label: "Grid",
        current_value: [
            g({
                type: "check",
                label: "Show density fields as a grid",
                current_value: config.grid?.by === "densityFields",
                func: (x) => {
                    config.grid = densityGridConfig(x);
                },
            }),
        ],
    });
}
