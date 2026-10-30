import { convertLegacyChartsSetting } from "@/lib/browserSettings";
import type DataStore from "@/datastore/DataStore";
import type { ChartType } from "./ChartTypes";
import { deckRadiusForClassicRadius } from "@/react/axisScaling";

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

/** The classic scatter's point size is in pixels, and the deck chart's is relative to the chart. */
const adaptClassicScatter = (config: any) => {
    adaptScatterBackgroundFilter(config);
    const size = Array.isArray(config.size) && config.size.length === 2 ? config.size : [300, 300];
    const sizing = deckRadiusForClassicRadius(config.radius ?? 2, size);
    if (sizing) Object.assign(config, sizing);
};

/**
 * Where a legacy chart type has a current equivalent that reads the same config keys.
 * Only these are converted: the Viv image viewer's config has no direct SpatialData.js equivalent.
 */
const legacyConversions: Record<string, LegacyConversion> = {
    wgl_scatter_plot: { to: "DeckContourScatter", adapt: adaptClassicScatter },
    // same chart class; `DeckContourScatter` only adds density settings, which are off by default
    wgl_scatter_plot_dev: { to: "DeckContourScatter", adapt: adaptScatterBackgroundFilter },
    // the sort direction is stored as `ascending` rather than `sortAsc`; the rest is shared
    table_chart: {
        to: "table_chart_react",
        adapt: (config) => {
            if (config.sort) config.sort = { columnId: config.sort.columnId, ascending: !!config.sort.sortAsc };
        },
    },
    // `param[2]` (the category column) becomes `contourParameter` when the deck chart loads
    // the Viv image viewer that is not a React chart; one whose region is SpatialData would want the
    // SpatialData.js viewer instead, which is not something this converts to
    viv_scatter_plot: {
        to: "VivMdvRegionReact",
        canConvert: (config, dataStore) => {
            const region = dataStore?.regions?.all_regions?.[config.region];
            return !!region?.viv_image && !region.spatial;
        },
    },
    density_scatter_plot: {
        to: "DeckContourScatter",
        canConvert: (config) => !config.viv,
        adapt: adaptClassicScatter,
    },
};

const findLegacyConversion = (config: { type: string }, dataStore?: DataStore) => {
    const conversion = legacyConversions[config.type];
    return conversion && conversion.canConvert?.(config, dataStore) !== false ? conversion : undefined;
};

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
