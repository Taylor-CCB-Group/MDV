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
};

const scatterBackgroundFilterIsUnset = (config: any) => {
    const column = config.background_filter?.column;
    return !column || column === "__none__";
};

/**
 * Where a legacy chart type has a current equivalent that reads the same config keys.
 * Only these are converted: the Viv image viewer's config has no direct SpatialData.js equivalent.
 */
const legacyConversions: Record<string, LegacyConversion> = {
    // the deck scatter ignores the classic charts' background filter, so a chart using one would change
    wgl_scatter_plot: { to: "DeckContourScatter", canConvert: scatterBackgroundFilterIsUnset },
    // same chart class; `DeckContourScatter` only adds density settings, which are off by default
    wgl_scatter_plot_dev: { to: "DeckContourScatter" },
    // `param[2]` (the category column) becomes `contourParameter` when the deck chart loads
    density_scatter_plot: {
        to: "DeckContourScatter",
        canConvert: (config) => !config.viv && scatterBackgroundFilterIsUnset(config),
    },
};

/**
 * If the "Convert legacy charts" browser setting is on and `config` is of a legacy type with a
 * current equivalent, changes its `type` in place so it loads as the equivalent.
 * @returns the type the chart was converted from, or undefined if it was left alone.
 */
export function convertLegacyChartConfig(config: { type: string }): string | undefined {
    if (convertLegacyChartsSetting.get() !== "on") return undefined;
    const conversion = legacyConversions[config.type];
    if (!conversion || conversion.canConvert?.(config) === false) return undefined;
    const from = config.type;
    config.type = conversion.to;
    return from;
}
