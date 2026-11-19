import type { ChartGridConfig } from "@/charts/schemas/ChartConfigSchema";
import type { DualContourLegacyConfig } from "../contour_state";
import { getChartArrayViewId } from "./chartArrayUtils";

/** Chart types that can show their density fields as a grid. */
const DENSITY_GRID_CHART_TYPES = new Set(["DeckContourScatter", "DeckDensity", "SpatialDataMdvRegionReact"]);

export function supportsDensityGridMode(chartType: string | undefined) {
    return typeof chartType === "string" && DENSITY_GRID_CHART_TYPES.has(chartType);
}

export type DensityGridModeConfig = Pick<DualContourLegacyConfig, "grid"> & {
    type?: string;
    dimension?: "2d" | "3d";
};

/** Whether the chart draws its density fields as a grid of cells rather than over its points. */
export function isDensityGridShown(config: DensityGridModeConfig) {
    // a chart with no dimension is 2D
    return supportsDensityGridMode(config.type) && config.dimension !== "3d" && config.grid?.by === "densityFields";
}

/**
 * The `grid` config that shows one cell per density field, or null for a single plot. The
 * current grid's other options are kept.
 */
export function densityGridConfig(shown: boolean, current?: ChartGridConfig | null): ChartGridConfig | null {
    return shown ? { ...current, by: "densityFields" } : null;
}

/** Whether the grid marks the pointer in every cell: on unless turned off. */
export function isGridCrosshairShown(grid: ChartGridConfig | null | undefined) {
    return grid?.crosshair ?? true;
}

/**
 * The `grid` config with its crosshair on or off. A new object: setting a key the config
 * object did not have would not be seen by MobX.
 */
export function withGridCrosshair(grid: ChartGridConfig | null | undefined, crosshair: boolean) {
    return grid ? { ...grid, crosshair } : grid;
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
