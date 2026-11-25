import { type RenderStack, walkRenderStackEntries } from "@spatialdata/layers";
import type { SpatialDataGridConfig } from "@/charts/schemas/ChartConfigSchema";
import { g } from "@/lib/utils";
import { renderStackEntryDisplayName } from "@/react/spatialdata/render_stack_display";
import { spatialDataUrlForRegion } from "@/react/spatialdata/shared_spatial_data";
import type { ChartArrayCell } from "../hooks/useChartArrayGrid";
import { COLOR_GRID_CHOICES, colorGridChoice, colorGridConfig, gridChoiceValue } from "./colorGridUtils";
import { layerFillGridChoice, layerFillGridChoices, layerFillGridConfig } from "./layerFillGridUtils";

/** A cell of a SpatialData chart drawn with a cell per region: the region and where its store is. */
export type SpatialRegionCell = ChartArrayCell & {
    /** The SpatialData object, as the converter named it (`region.spatial.file`). */
    file: string;
    url: string;
    coordinateSystem: string;
};

type RegionGridModeConfig = {
    type?: string;
    grid?: { by: string } | null;
};

/** Whether the SpatialData chart draws a cell per region rather than its one region. */
export function isRegionGridShown(config: RegionGridModeConfig) {
    return config.type === "SpatialDataMdvRegionReact" && config.grid?.by === "region";
}

function spatialMetadata(region: unknown): { file: string; coordinateSystem: string } | null {
    if (!region || typeof region !== "object" || !("spatial" in region)) return null;
    const { spatial } = region;
    if (!spatial || typeof spatial !== "object") return null;
    const file = "file" in spatial ? spatial.file : undefined;
    const coordinateSystem = "coordinate_system" in spatial ? spatial.coordinate_system : undefined;
    if (typeof file !== "string" || !file || typeof coordinateSystem !== "string" || !coordinateSystem) return null;
    return { file, coordinateSystem };
}

/** A cell for each region with SpatialData, in the datasource's order. */
export function spatialRegionCells(allRegions: Record<string, unknown> | undefined): SpatialRegionCell[] {
    if (!allRegions) return [];
    return Object.entries(allRegions).flatMap(([key, region]) => {
        const metadata = spatialMetadata(region);
        const url = spatialDataUrlForRegion(region);
        if (!metadata || !url) return [];
        return [{ key, label: key, url, ...metadata }];
    });
}

/** A region cell's view: where it is centred, in the region's coordinate system, and its zoom. */
export type RegionViewState = { target: [number, number]; zoom: number };

/** A saved view state as a region cell's, with what is missing filled in. */
export function toRegionViewState(state: { target?: ArrayLike<number>; zoom?: unknown }): RegionViewState {
    const target = state.target ?? [];
    return { target: [target[0] ?? 0, target[1] ?? 0], zoom: typeof state.zoom === "number" ? state.zoom : 0 };
}

/**
 * Every region's view after one cell's moved to `next`. Linked, the other cells move as it did on
 * screen: the same change of zoom, and the same pan in pixels, so each keeps its own framing and a
 * zoom about the pointer zooms each about the same point of its cell. Cells without a view yet are
 * left to be fitted.
 */
export function moveRegionViews(
    states: Readonly<Record<string, RegionViewState>>,
    region: string,
    next: RegionViewState,
    linked: boolean,
): Record<string, RegionViewState> {
    const previous = states[region];
    const moved: Record<string, RegionViewState> = { ...states, [region]: next };
    if (!linked || !previous) return moved;
    const zoomChange = next.zoom - previous.zoom;
    // the pan in pixels at the new zoom, the same in every cell
    const panX = (next.target[0] - previous.target[0]) * 2 ** next.zoom;
    const panY = (next.target[1] - previous.target[1]) * 2 ** next.zoom;
    for (const [key, state] of Object.entries(states)) {
        if (key === region) continue;
        const zoom = state.zoom + zoomChange;
        const scale = 2 ** zoom;
        moved[key] = { target: [state.target[0] + panX / scale, state.target[1] + panY / scale], zoom };
    }
    return moved;
}

/**
 * The SpatialData chart's grid settings: one view, a cell per density field, per colour column,
 * per category of the first colour column, per fill colour column (or category) of a shapes or
 * labels layer, or per region.
 */
export function getSpatialDataGridSettingsFolder(config: {
    grid?: SpatialDataGridConfig | null;
    renderStack?: RenderStack;
}) {
    const choices = [
        { name: "One view", value: "none" },
        { name: "A cell per density field", value: "densityFields" },
        ...COLOR_GRID_CHOICES,
        ...layerFillGridChoices(fillColorLayers(config.renderStack)),
        { name: "A cell per region", value: "region" },
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
                    const layerGrid = layerFillGridChoice(value);
                    if (colorSplit) config.grid = colorGridConfig(colorSplit);
                    else if (layerGrid) config.grid = layerFillGridConfig(layerGrid.layer, layerGrid.split);
                    else if (value === "densityFields") config.grid = { by: "densityFields" };
                    else if (value === "region")
                        config.grid = {
                            by: "region",
                            hide_empty: config.grid?.hide_empty,
                            link_views: config.grid?.link_views,
                        };
                    else config.grid = null;
                },
            }),
            g({
                type: "check",
                label: "Leave out regions with no rows after filtering",
                current_value: config.grid?.hide_empty ?? false,
                func: (x: boolean) => {
                    if (config.grid?.by === "region") config.grid = { ...config.grid, hide_empty: x || undefined };
                },
            }),
            g({
                type: "check",
                label: "Pan and zoom every region together",
                current_value: config.grid?.link_views ?? false,
                func: (x: boolean) => {
                    if (config.grid?.by === "region") config.grid = { ...config.grid, link_views: x || undefined };
                },
            }),
        ],
    });
}

/** The stack's shapes and labels layers, which can be coloured by a column: their ids and names. */
export function fillColorLayers(stack: RenderStack | undefined): { id: string; name: string }[] {
    if (!stack) return [];
    return Array.from(walkRenderStackEntries(stack.entries)).flatMap((entry) =>
        entry.kind === "spatial" && (entry.source.elementType === "shapes" || entry.source.elementType === "labels")
            ? [{ id: entry.id, name: renderStackEntryDisplayName(entry, stack) }]
            : [],
    );
}
