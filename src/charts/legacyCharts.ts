import type DataStore from "@/datastore/DataStore";
import { convertLegacyChartsSetting } from "@/lib/browserSettings";
import { deckRadiusForClassicRadius, deckRadiusForPixelRadius, radiusReferenceRange3D } from "@/react/axisScaling";
import { adaptMultiLineConfig } from "@/react/components/multiline/multiLineConfig";
import { fitOrbitView, orbitViewFromLegacyCamera, validAxisScales } from "@/react/orbitView";
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
    adapt?: (config: any, dataStore?: DataStore) => void;
    /**
     * Convert charts of this type whenever a view loads, not only with the "Convert legacy charts"
     * browser setting on: for a legacy chart that its replacement now stands in for everywhere.
     */
    always?: boolean;
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

/**
 * The classic chart's title bar (31px) and the margins round its plot (20px): what `size` has that
 * the plot itself does not, as measured on the classic chart.
 */
const CLASSIC_CHROME_HEIGHT = 50;

/** The `[min, max]` of each of a chart's x, y and z columns, if the data store knows them. */
function columnRanges(config: any, dataStore?: DataStore) {
    const param: unknown[] = Array.isArray(config.param) ? config.param.slice(0, 3) : [];
    if (!dataStore || param.length < 3) return undefined;
    try {
        const ranges = param.map((p) => (typeof p === "string" ? dataStore.getMinMaxForColumn(p) : undefined));
        if (ranges.every((r): r is [number, number] => Array.isArray(r) && r.every(Number.isFinite))) return ranges;
    } catch {
        // a column with neither min/max metadata nor loaded data: leave the point size as it was
    }
    return undefined;
}

/**
 * The classic (WebGL) 3D scatter as the deck chart. Its camera becomes an equivalent `viewState`,
 * with y up as the classic chart had it (see `orbitViewFromLegacyCamera`), and `axis_scales`
 * carries over as it is. Point size was a pixel diameter that stayed the same at any zoom; it
 * becomes the deck size that draws points that size at the saved view, where the data store has the
 * column ranges to work that out.
 *
 * A `viewState` from before the classic chart came back (when this type was briefly the deck chart)
 * is dropped: since then the classic chart has drawn it, and saved its own `camera`.
 */
const adaptClassic3DScatter = (config: any, dataStore?: DataStore) => {
    const scales = validAxisScales(config.axis_scales);
    const [width, height] = chartSize(config);
    const plotWidth = Math.max(width - 20, 1);
    const plotHeight = Math.max(height - CLASSIC_CHROME_HEIGHT, 1);
    const camera = config.camera;
    const c = config.center;
    const view =
        camera && Array.isArray(c) && c.length === 3
            ? orbitViewFromLegacyCamera(camera, [c[0], c[1], c[2]], plotHeight)
            : undefined;
    // the classic 3D chart never applied a background or category filter (only the classic 2D one
    // did), so one carried over in its config would start hiding rows the chart showed
    for (const key of ["camera", "center", "brush", "viewState", "background_filter", "category_filters"]) {
        delete config[key];
    }
    config.dimension = "3d";
    config.orbit_axis = "Y";
    config.axis_scales = scales;
    if (view) config.viewState = view;

    const ranges = columnRanges(config, dataStore);
    if (!ranges) return;
    // where the chart starts, which is what the classic size is matched at
    const zoom =
        view?.zoom ??
        fitOrbitView(
            {
                min: ranges.map(([lo], i) => lo * scales[i]),
                max: ranges.map(([, hi], i) => hi * scales[i]),
            },
            plotWidth,
            plotHeight,
        ).zoom;
    const referenceRange = radiusReferenceRange3D(ranges.map(([lo, hi]) => hi - lo));
    if (!(referenceRange > 0)) return;
    const sizing = deckRadiusForPixelRadius((config.radius ?? 5) / 2, 2 ** zoom, referenceRange);
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
    // a pixel-sized point and a camera of its own, which `adaptClassic3DScatter` converts
    wgl_3d_scatter_plot: [{ to: "DeckScatter3D", adapt: adaptClassic3DScatter, always: true }],
    // the same chart class under an older name
    wgl_3d_scatter_plot_dev: [{ to: "DeckScatter3D", adapt: (config) => void (config.dimension = "3d"), always: true }],
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
    // `param: [value, category]` becomes `values` and `categories`; the rest reads the same
    multi_line_chart: [
        { to: "multi_line_chart_react", adapt: (config) => Object.assign(config, adaptMultiLineConfig(config)) },
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
    conversion.adapt?.(config, dataStore);
    return from;
}

/**
 * Like {@link convertLegacyChartConfigNow}, but only if the "Convert legacy charts" browser setting is on,
 * or the chart's type is one that always converts: used when charts are loaded.
 */
export function convertLegacyChartConfig(config: { type: string }, dataStore?: DataStore): string | undefined {
    if (convertLegacyChartsSetting.get() !== "on" && !findLegacyConversion(config, dataStore)?.always) {
        return undefined;
    }
    return convertLegacyChartConfigNow(config, dataStore);
}
