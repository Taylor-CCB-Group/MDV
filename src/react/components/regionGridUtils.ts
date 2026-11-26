import { type RenderStack, walkRenderStackEntries } from "@spatialdata/layers";
import type { SpatialDataGridConfig } from "@/charts/schemas/ChartConfigSchema";
import type { FieldSpecs } from "@/lib/columnTypeHelpers";
import { g } from "@/lib/utils";
import { renderStackEntryDisplayName } from "@/react/spatialdata/render_stack_display";
import { spatialDataUrlForRegion } from "@/react/spatialdata/shared_spatial_data";
import { availableColorGridSplit, type ColumnTypes } from "../colorFields";
import type { ChartArrayCell } from "../hooks/useChartArrayGrid";
import {
    type ColorGridSplit,
    colorGridChoice,
    colorGridChoices,
    colorGridConfig,
    colorGridSplit,
    type GridChoice,
    gridChoiceValue,
    ONE_VIEW,
} from "./colorGridUtils";
import { fillColumnName, layerFillGridChoice, layerFillGridChoices, layerFillGridConfig } from "./layerFillGridUtils";
import { getGreyOtherCategoriesSetting } from "./scatterGridSettings";

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

const SPLITS: readonly ColorGridSplit[] = ["columns", "categories"];

type SpatialDataGridSettingsConfig = {
    type?: string;
    grid?: SpatialDataGridConfig | null;
    renderStack?: RenderStack;
    color_fields?: FieldSpecs;
    grey_other_categories?: boolean | null;
};

/**
 * Each layer's colour grid choices: the split its colour column allows, and the chart's current
 * grid, so the dropdown can show it.
 */
function layerGridChoices(config: SpatialDataGridSettingsConfig, layerGridSplit: LayerGridSplit | undefined) {
    const { grid } = config;
    const current =
        grid?.by === "layer_fill" && grid.layer
            ? { layer: grid.layer, split: grid.split === "categories" ? "categories" : "columns" }
            : undefined;
    return fillColorLayers(config.renderStack).flatMap((layer) => {
        const available = layerGridSplit?.(layer);
        return SPLITS.filter(
            (split) =>
                !layerGridSplit || split === available || (current?.layer === layer.id && current.split === split),
        ).map((split) => ({ id: layer.id, name: layer.name, split }));
    });
}

/**
 * How a grid of a layer's fill colour would be split, by the type of its column in the table
 * annotating its element, or null when it cannot be one (no table, or no column of a type to
 * split). Absent when that cannot be told yet: every layer is offered either way.
 */
type LayerGridSplit = (layer: FillColorLayer) => ColorGridSplit | null;

/**
 * The grids the SpatialData chart can show, besides one view: a cell per density field, a grid of
 * the chart's colour or of a shapes or labels layer's fill colour (a cell per colour column, or per
 * category of a categorical one), or a cell per region.
 */
export function spatialDataGridChoices(
    config: SpatialDataGridSettingsConfig,
    columns: ColumnTypes,
    layerGridSplit?: LayerGridSplit,
): GridChoice[] {
    return [
        { name: "A cell per density field", value: "densityFields" },
        ...colorGridChoices(availableColorGridSplit(config, columns), colorGridSplit(config)),
        ...layerFillGridChoices(layerGridChoices(config, layerGridSplit)),
        { name: "A cell per region", value: "region" },
    ];
}

/** Shows the grid a choice names (`spatialDataGridChoices`), or one view for "none". */
export function setSpatialDataGrid(config: SpatialDataGridSettingsConfig, value: string) {
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
}

/** The SpatialData chart's grid settings: one view, or a grid (`spatialDataGridChoices`). */
export function getSpatialDataGridSettingsFolder(
    config: SpatialDataGridSettingsConfig,
    columns: ColumnTypes,
    layerGridSplit?: LayerGridSplit,
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
                // the colour choices follow the colour columns' types
                get values(): [GridChoice[], string, string] {
                    return [[ONE_VIEW, ...spatialDataGridChoices(config, columns, layerGridSplit)], "name", "value"];
                },
                set values(_values) {
                    // made from the config
                },
                func: (value: string) => setSpatialDataGrid(config, value),
            }),
            getGreyOtherCategoriesSetting(config),
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

/** A shapes or labels layer of the stack, which can be coloured by a column if a table annotates it. */
export type FillColorLayer = {
    id: string;
    name: string;
    elementType: "shapes" | "labels";
    elementKey: string;
    /** The column it is coloured by (the first, for several), if any. */
    fillColumn?: string;
};

/** The stack's shapes and labels layers: their ids, names and elements. */
export function fillColorLayers(stack: RenderStack | undefined): FillColorLayer[] {
    if (!stack) return [];
    return Array.from(walkRenderStackEntries(stack.entries)).flatMap((entry) => {
        if (entry.kind !== "spatial") return [];
        const { elementType, elementKey } = entry.source;
        if (elementType !== "shapes" && elementType !== "labels") return [];
        const name = renderStackEntryDisplayName(entry, stack);
        return [{ id: entry.id, name, elementType, elementKey, fillColumn: fillColumnName(entry.props) }];
    });
}
