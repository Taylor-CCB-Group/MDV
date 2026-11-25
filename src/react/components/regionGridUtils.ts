import type { SpatialDataGridConfig } from "@/charts/schemas/ChartConfigSchema";
import { g } from "@/lib/utils";
import { spatialDataUrlForRegion } from "@/react/spatialdata/shared_spatial_data";
import type { ChartArrayCell } from "../hooks/useChartArrayGrid";

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

/** A region's rows, and where each is in the array they were taken from. */
export type RegionRows = { rows: Uint32Array; positions: Uint32Array };

/**
 * Splits rows by region, keeping their order. `codes` is the region column's data (an index into
 * `values` per row). Regions with no rows are left out.
 */
export function partitionRowsByRegion(
    rows: ArrayLike<number>,
    codes: ArrayLike<number>,
    values: readonly string[],
): Map<string, RegionRows> {
    const counts = new Uint32Array(values.length);
    for (let i = 0; i < rows.length; i++) {
        const code = codes[rows[i]];
        if (code < values.length) counts[code]++;
    }
    const byCode = Array.from(counts, (count) =>
        count > 0 ? { rows: new Uint32Array(count), positions: new Uint32Array(count) } : null,
    );
    const filled = new Uint32Array(values.length);
    for (let i = 0; i < rows.length; i++) {
        const row = rows[i];
        const code = codes[row];
        const target = byCode[code];
        if (!target) continue;
        const at = filled[code]++;
        target.rows[at] = row;
        target.positions[at] = i;
    }
    const result = new Map<string, RegionRows>();
    byCode.forEach((target, code) => {
        if (target) result.set(values[code], target);
    });
    return result;
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

type GridChoice = "none" | SpatialDataGridConfig["by"];

/**
 * The SpatialData chart's grid settings: one view, a cell per density field, or a cell per region.
 */
export function getSpatialDataGridSettingsFolder(config: { grid?: SpatialDataGridConfig | null }) {
    const choices: { name: string; value: GridChoice }[] = [
        { name: "One view", value: "none" },
        { name: "A cell per density field", value: "densityFields" },
        { name: "A cell per region", value: "region" },
    ];
    return g({
        type: "folder",
        label: "Grid",
        current_value: [
            g({
                type: "dropdown",
                label: "Show",
                current_value: config.grid?.by ?? "none",
                values: [choices, "name", "value"],
                func: (value: string) => {
                    if (value === "densityFields") config.grid = { by: "densityFields" };
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
