import type DataStore from "@/datastore/DataStore";
import { convertLegacyChartsSetting } from "@/lib/browserSettings";
import { deckRadiusForClassicRadius } from "@/react/axisScaling";
import type { ChartType } from "./ChartTypes";

/** Whether a chart type is hidden from the Add Chart dialog as superseded, for this data store. */
export function isLegacyChartType(type: Pick<ChartType<any>, "legacy">, dataStore: DataStore) {
    return typeof type.legacy === "function" ? type.legacy(dataStore) : !!type.legacy;
}

type LegacyConversion = {
    /** The chart type that replaces it. */
    to: string;
    /** False where the replacement would not show what the saved chart does. */
    canConvert?: (config: any, dataStore?: DataStore) => boolean;
    /** Rewrites config keys that the replacement reads differently. */
    adapt?: (config: any) => void;
};

/**
 * The classic scatter keeps the selected values in `categories` (older charts in `category`) and
 * counts a column of `__none__` or no categories as unset; the deck scatter reads `category`.
 */
const adaptScatterBackgroundFilter = (config: any) => {
    const filter = config.background_filter;
    if (!filter) return;
    const categories: unknown[] = filter.categories?.length ? filter.categories : [filter.category];
    const category = categories.filter((x) => typeof x === "string" && x !== "");
    if (!filter.column || filter.column === "__none__" || category.length === 0) {
        config.background_filter = undefined;
    } else {
        config.background_filter = { column: filter.column, category };
    }
};

/** The chart's size, from `size` or else from `location`, which `addChart` turns into a size and position. */
const chartSize = (config: any): [number, number] => {
    if (Array.isArray(config.size) && config.size.length === 2) return config.size as [number, number];
    const l = config.location;
    if (l && typeof l.width === "number" && typeof l.height === "number") {
        const b = 5;
        return [l.width * 90 + l.width * b - b, l.height * 40 + l.height * b - b];
    }
    return [300, 300];
};

/** The classic scatter's point size is in pixels, and the deck chart's is relative to the chart. */
const adaptClassicScatter = (config: any) => {
    adaptScatterBackgroundFilter(config);
    const size = chartSize(config);
    const sizing = deckRadiusForClassicRadius(config.radius ?? 2, size);
    if (sizing) Object.assign(config, sizing);
};

/** The region a Viv chart shows, as `dataStore.regions` describes it. */
const vivRegion = (config: any, dataStore?: DataStore) => dataStore?.regions?.all_regions?.[config.region];

/** A Viv image region that is also SpatialData, which the SpatialData.js viewer can show. */
const isSpatialDataVivRegion = (config: any, dataStore?: DataStore) => {
    const region = vivRegion(config, dataStore);
    return !!region?.viv_image && !!region.spatial;
};

/**
 * The Viv chart's saved channel settings as the SpatialData.js image layer's: the same per-channel
 * arrays, with tone in `vivLayerProps`. Undefined if it saved none, so the layer starts with defaults.
 */
const vivChannelsToImageLayer = (channelsStore: any) => {
    if (!channelsStore || typeof channelsStore !== "object") return undefined;
    const channels: Record<string, unknown> = {};
    for (const key of ["colors", "contrastLimits", "channelsVisible", "selections"]) {
        if (Array.isArray(channelsStore[key])) channels[key] = channelsStore[key];
    }
    const vivLayerProps: Record<string, unknown> = {};
    for (const key of ["brightness", "contrast"]) {
        if (Array.isArray(channelsStore[key]) && channelsStore[key].length) vivLayerProps[key] = channelsStore[key];
    }
    if (!Object.keys(channels).length && !Object.keys(vivLayerProps).length) return undefined;
    return {
        ...(Object.keys(channels).length ? { channels } : {}),
        ...(Object.keys(vivLayerProps).length ? { vivLayerProps } : {}),
    };
};

/**
 * The SpatialData.js viewer reads the same scatter, region and filter settings as the Viv viewer, but
 * draws images from a render stack rather than `viv`. Left without one, it adds the store's first image
 * on load, with the Viv channel settings from `defaultImageLayer`. The view position moves to the
 * top-level `viewState`: both viewers draw the scatter in the same frame.
 */
const adaptVivToSpatialData = (config: any) => {
    const vivView = config.viv?.viewerStore?.viewState;
    const defaultImageLayer = vivChannelsToImageLayer(config.viv?.channelsStore);
    for (const key of ["viv", "roi", "json", "showJson", "background_image", "renderStack", "viewState"]) {
        delete config[key];
    }
    if (Array.isArray(vivView?.target)) {
        config.viewState = { target: [...vivView.target], zoom: vivView.zoom };
    }
    if (defaultImageLayer) config.defaultImageLayer = defaultImageLayer;
};

/**
 * Where a legacy chart type has current equivalents, the first that works for a chart's config is
 * the one it converts to.
 */
const legacyConversions: Record<string, LegacyConversion[]> = {
    wgl_scatter_plot: [{ to: "DeckContourScatter", adapt: adaptClassicScatter }],
    // same chart class; `DeckContourScatter` only adds density settings, which are off by default
    wgl_scatter_plot_dev: [{ to: "DeckContourScatter", adapt: adaptScatterBackgroundFilter }],
    // the sort direction is stored as `ascending` rather than `sortAsc`; the rest is shared
    table_chart: [
        {
            to: "table_chart_react",
            adapt: (config) => {
                if (config.sort) config.sort = { columnId: config.sort.columnId, ascending: !!config.sort.sortAsc };
            },
        },
    ],
    // `param[2]` (the category column) becomes `contourParameter` when the deck chart loads
    viv_scatter_plot: [
        { to: "SpatialDataMdvRegionReact", canConvert: isSpatialDataVivRegion, adapt: adaptVivToSpatialData },
        {
            to: "VivMdvRegionReact",
            canConvert: (config, dataStore) => !!vivRegion(config, dataStore)?.viv_image,
        },
    ],
    // still the viewer for a Viv image region that is not SpatialData
    VivMdvRegionReact: [
        { to: "SpatialDataMdvRegionReact", canConvert: isSpatialDataVivRegion, adapt: adaptVivToSpatialData },
    ],
    density_scatter_plot: [
        {
            to: "DeckContourScatter",
            canConvert: (config) => !config.viv,
            adapt: adaptClassicScatter,
        },
    ],
};

const findLegacyConversion = (config: { type: string }, dataStore?: DataStore) =>
    legacyConversions[config.type]?.find((conversion) => conversion.canConvert?.(config, dataStore) !== false);

/** The chart type that `config` could be converted to, if its type is legacy and that would work. */
export function getLegacyConversionTarget(config: { type: string }, dataStore?: DataStore): string | undefined {
    return findLegacyConversion(config, dataStore)?.to;
}

/**
 * Changes `config`'s `type` in place to its current equivalent, if there is one that works for it.
 * @returns the type the chart was converted from, or undefined if it was left alone.
 */
export function convertLegacyChartConfigNow(config: { type: string }, dataStore?: DataStore): string | undefined {
    const conversion = findLegacyConversion(config, dataStore);
    if (!conversion) return undefined;
    const from = config.type;
    config.type = conversion.to;
    conversion.adapt?.(config);
    return from;
}

/**
 * Like {@link convertLegacyChartConfigNow}, but only if the "Convert legacy charts" browser setting is on:
 * used when charts are loaded.
 */
export function convertLegacyChartConfig(config: { type: string }, dataStore?: DataStore): string | undefined {
    if (convertLegacyChartsSetting.get() !== "on") return undefined;
    return convertLegacyChartConfigNow(config, dataStore);
}
