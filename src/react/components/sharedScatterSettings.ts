import type BaseChart from "@/charts/BaseChart";
import getTooltipSettings from "@/charts/dialogs/utils/TooltipSettingsGui";
import { g } from "@/lib/utils";
import { getDensitySettings } from "../contour_state";
import { scatterDefaults, type ScatterPlotConfig } from "../scatter_state";

type SharedScatterSettingsOptions<C extends ScatterPlotConfig> = {
    chart?: BaseChart<C>;
    includeDensitySettings?: boolean;
    includePointShape?: boolean;
    includeTooltip?: boolean;
    includeZoomOnFilter?: boolean;
};

export function getSharedScatterSettings<C extends ScatterPlotConfig>(
    config: C,
    {
        chart,
        includeDensitySettings = false,
        includePointShape = false,
        includeTooltip = true,
        includeZoomOnFilter = true,
    }: SharedScatterSettingsOptions<C> = {},
) {
    const settings = [];

    if (includeTooltip) {
        settings.push(getTooltipSettings(config));
    }

    if (includePointShape) {
        settings.push(
            g({
                type: "dropdown",
                label: "Shape",
                current_value: config.point_shape,
                values: [["circle", "square", "gaussian"]],
                func: (value) => {
                    // @ts-ignore existing GuiSpec dropdowns are stringly typed
                    config.point_shape = value;
                },
            }),
        );
    }

    settings.push(
        g({
            type: "radiobuttons",
            label: "course radius",
            current_value: `${config.course_radius || 1}`,
            choices: [
                ["0.1", "0.1"],
                ["1", "1"],
                ["10", "10"],
                ["100", "100"],
            ],
            func: (value) => {
                config.course_radius = Number.parseFloat(value);
            },
        }),
        g({
            type: "slider",
            label: "radius",
            current_value: config.radius || 5,
            min: 0,
            max: 20,
            continuous: true,
            func: (value) => {
                config.radius = value;
            },
        }),
        g({
            type: "slider",
            label: "opacity",
            current_value: Math.sqrt(config.opacity || scatterDefaults.opacity),
            min: 0,
            max: 1,
            continuous: true,
            func: (value) => {
                config.opacity = value * value;
            },
        }),
    );

    if (includeZoomOnFilter) {
        settings.push(
            g({
                type: "check",
                label: "zoom on filter",
                current_value: config.zoom_on_filter || false,
                func: (value) => {
                    config.zoom_on_filter = value;
                },
            }),
        );
    }

    if (includeDensitySettings && chart) {
        settings.push(getDensitySettings(config));
    }

    return settings;
}
