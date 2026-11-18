import type { DualContourLegacyConfig } from "../contour_state";
import { getChartArrayViewId } from "./chartArrayUtils";

/** Chart types that can show their density fields as a grid. */
const DENSITY_GRID_CHART_TYPES = new Set(["DeckContourScatter", "DeckDensity"]);

export function supportsDensityGridMode(chartType: string | undefined) {
    return typeof chartType === "string" && DENSITY_GRID_CHART_TYPES.has(chartType);
}

export type DensityGridModeConfig = Pick<DualContourLegacyConfig, "density_mode"> & {
    type?: string;
    dimension?: "2d" | "3d";
};

/** Whether the chart draws its density fields as a grid of cells rather than over its points. */
export function isDensityGridShown(config: DensityGridModeConfig) {
    return supportsDensityGridMode(config.type) && config.dimension === "2d" && config.density_mode === "grid";
}

const DENSITY_GRID_PREFIX = "density-grid";

/** The view id of a grid cell. Layers that belong to one cell take it as their id's prefix. */
export function getDensityGridViewId(chartId: string, fieldId: string, index: number) {
    return getChartArrayViewId(chartId, fieldId, index, DENSITY_GRID_PREFIX);
}

/**
 * The grid's `layerFilter`. A layer whose id starts with a cell's view id (as do its
 * sublayers', `<parent>-<child>`) draws only in that cell; every other layer draws in all
 * of them.
 */
export function shouldDrawInDensityGridCell(layerId: string, viewId: string) {
    if (!layerId.startsWith(`${DENSITY_GRID_PREFIX}-`)) return true;
    return layerId === viewId || layerId.startsWith(`${viewId}-`);
}
