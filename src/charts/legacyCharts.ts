import { convertLegacyChartsSetting } from "@/lib/browserSettings";
import type DataStore from "@/datastore/DataStore";
import type { ChartType } from "./ChartTypes";

/** Whether a chart type is hidden from the Add Chart dialog as superseded, for this data store. */
export function isLegacyChartType(type: Pick<ChartType<any>, "legacy">, dataStore: DataStore) {
    return typeof type.legacy === "function" ? type.legacy(dataStore) : !!type.legacy;
}

type LegacyConversion = {
    /** The chart type that replaces it. */
    to: string;
    /** False where the replacement would not show what the saved chart does. */
    canConvert?: (config: any) => boolean;
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

/**
 * Where a legacy chart type has a current equivalent that reads the same config keys.
 * Only these are converted: the Viv image viewer's config has no direct SpatialData.js equivalent.
 */
const legacyConversions: Record<string, LegacyConversion> = {
    wgl_scatter_plot: { to: "DeckContourScatter", adapt: adaptScatterBackgroundFilter },
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
    density_scatter_plot: {
        to: "DeckContourScatter",
        canConvert: (config) => !config.viv,
        adapt: adaptScatterBackgroundFilter,
    },
};

/**
 * Changes `config`'s `type` in place to its current equivalent, if there is one that works for it.
 * @returns the type the chart was converted from, or undefined if it was left alone.
 */
export function convertLegacyChartConfigNow(config: { type: string }): string | undefined {
    const conversion = findLegacyConversion(config);
    if (!conversion) return undefined;
const findLegacyConversion = (config: { type: string }) => {
    const conversion = legacyConversions[config.type];
    return conversion && conversion.canConvert?.(config) !== false ? conversion : undefined;
};

/** The chart type that `config` could be converted to, if its type is legacy and that would work. */
export function getLegacyConversionTarget(config: { type: string }): string | undefined {
    return findLegacyConversion(config)?.to;
}

    const from = config.type;
    config.type = conversion.to;
    conversion.adapt?.(config);
    return from;
}

/**
 * Like {@link convertLegacyChartConfigNow}, but only if the "Convert legacy charts" browser setting is on:
 * used when charts are loaded.
 */
export function convertLegacyChartConfig(config: { type: string }): string | undefined {
    if (convertLegacyChartsSetting.get() !== "on") return undefined;
    return convertLegacyChartConfigNow(config);
}
