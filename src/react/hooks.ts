import { useEffect, useLayoutEffect, useMemo, useState } from "react";
import { loadOmeTiff, getChannelStats } from "@hms-dbmi/viv";
import { useChart, useDataStore } from "./context";
import type { OME_TIFF } from "./components/avivatorish/state";
import { getProjectURL, loadColumn } from "../dataloaders/DataLoaderUtil";
import { getRandomString } from "../utilities/Utilities";
import { action } from "mobx";
import type { CategoricalDataType, DataColumn } from "../charts/charts";
import type { VivRoiConfig } from "./components/VivMDVReact";

/**
 * Get the chart's config.
 * 
 * Must be used within a ChartContext.
 * 
 * Provided type parameter is not checked - in future it could probably
 * be inferred from the chart type.
 */
export function useConfig<T>() {
    const { config } = useChart();
    return config as T; //todo: strict/inferred typing
}

export function useChartSize() {
    const chart = useChart();
    // return chart.config.size; // not so well behaved?
    const div = chart.contentDiv;
    const [size, setSize] = useState([div.clientWidth, div.clientHeight]);
    useLayoutEffect(() => {
        const resize = () => {
            setSize([div.clientWidth, div.clientHeight]);
        };
        const observer = new ResizeObserver(resize);
        observer.observe(div);
        return () => observer.unobserve(div);
    }, [div]);
    return size;
}

/**
 * Get the chart's ID.
 * 
 * Must be used within a ChartContext.
 */
export function useChartID(): string {
    const chart = useChart();
    if (!chart.config.id) {
        // we were hitting this because of the way BaseReactChart was still using original config
        // after super constructor had used a copy of it to set up the chart...
        console.assert(chart.config.id, 'chart.config.id should not be undefined');
        action(() => {
            chart.config.id = getRandomString();
        })();
    }
    return chart.config.id;
}

/** Get the document the chart is currently assigned to, which changes when popped-out.
 * Used for mouse events etc that were previously on `window` to allow dragging out of the chart itself.
 * @deprecated - use `useOuterContainer` instead.
 */
export function useChartDoc() {
    const chart = useChart();
    return chart.__doc__;
}

export function useParamColumns(): DataColumn<any>[] {
    const chart = useChart();
    const { columnIndex } = chart.dataStore;
    const columns = useMemo(() => {
        const param = chart.config.param;
        if (!param) return [];
        if (typeof chart.config.param === 'string') return [columnIndex[chart.config.param]];
        // we should make sure they are loaded as well...
        return chart.config.param.map(name => columnIndex[name])
    }, [chart.config.param, columnIndex]);
    return columns;
}

export function useNamedColumn(name: string): { column: DataColumn<any>, isLoaded: boolean } {
    const chart = useChart();
    const { columnIndex } = chart.dataStore;
    const [isLoaded, setIsLoaded] = useState(false);
    const column = columnIndex[name];
    useEffect(() => {
        loadColumn(chart.dataStore.name, name).then(() => setIsLoaded(true));
    }, [name, chart.dataStore]);
    return { column, isLoaded };
}

/** If the chart in current context has an associated region, referred to by the key `config.region`, this should return it
 * such that relevant information (like image URL, associated geojson layers...) can be retrieved.
 */
export function useRegion() {
    const config = useConfig<VivRoiConfig>();
    const { regions } = useDataStore();
    const region = useMemo(() => {
        if (!regions) return undefined;
        //we could add a dash of zod here?
        //at the moment, this just returns 'any'.
        return regions.all_regions[config.region];
    }, [config.region, regions]);
    return region;
}

/**
 * Get a {Uint32Array} of the currently filtered indices.
 * When the selection changes, this will asynchronously update.
 * All users of the same data store (on a per-chart basis) share a reference to the same array.
 * -- change properties/settings so that we can more dynamically select.
 */
export function useFilteredIndices() {
    // in the case of region data, it should be filtered by that as well...
    // I really want to sort out how I use types here...
    const config = useConfig<VivRoiConfig>();
    const filterColumn = config.background_filter?.column;
    const dataStore = useDataStore();
    const [filteredIndices, setFilteredIndices] = useState(new Uint32Array());
    const [filteredOutIndices, setFilteredOutIndices] = useState(new Uint32Array());
    // biome-ignore lint/correctness/useExhaustiveDependencies: shouldn't be ignoring this, but some deps don't tend to change as of now.
    useEffect(() => {
        // return
        let cancelled = false;
        if (!filterColumn) return;
        const indexPromise = dataStore.getFilteredIndices();
        //todo maybe make more use of deck.gl category filters once we update to new version
        const catFilters = [config.background_filter, ...config.category_filters.filter(f => f.category !== 'all')];
        const catColumns = catFilters.map(f => f.column);
        const colPromise = window.mdv.chartManager?._getColumnsAsync(dataStore.name, catColumns);
        Promise.all([indexPromise, colPromise]).then(([indices]) => {
            if (cancelled) return;
            if (filterColumn) {
                const cols = catFilters.map(({ column }) => dataStore.columnIndex[column]);
                const filterValue = config.background_filter?.category;
                if (filterValue) {
                    //const filterIndex = col.values.indexOf(filterValue);
                    const filterIndex = catFilters.map(f => {
                        if (Array.isArray(f.category)) return f.category.map(c => dataStore.columnIndex[f.column].values.indexOf(c)) as number[];
                        return dataStore.columnIndex[f.column].values.indexOf(f.category) as number;
                    });
                    try {
                        // const filteredIndices = indices.filter(i => col.data[i] === filterIndex);
                        const filteredIndices = indices.filter(i => catFilters.every((_, j) => {
                            const f = filterIndex[j];
                            if (typeof f === 'number') return f === cols[j].data[i];
                            return f.some(fi => cols[j].data[i] === fi);
                        }));
                        setFilteredIndices(filteredIndices);
                        // thinking about allowing gray-out of non-selected points... should be optional
                        // const filteredOutIndices = indices.filter(i => col.data[i] !== filterIndex);
                        // setFilteredOutIndices(filteredOutIndices);
                    } catch (e) {
                        console.error('error filtering indices', e);
                        return;
                    }
                    return;
                }
            }
            setFilteredIndices(indices);
        });
        // should I have a cleanup function to cancel the promise if it's not resolved
        // by the time the effect is triggered again?
        return () => {
            // if (!finished) console.log('filtered indices promise cancelled');
            cancelled = true;
        }

        // using _filteredIndicesPromise as a dependency is working reasonably well,
        // but possibly needs a bit more thought.
    }, [dataStore._filteredIndicesPromise, filterColumn, config.background_filter, config.category_filters]);
    return filteredIndices;
}

export function useCategoryFilterIndices(contourParameter: DataColumn<CategoricalDataType>, category: string | string[]) {
    const data = useFilteredIndices();
    //todo handle multitext / tags properly.
    const categoryValueIndex = useMemo(() => {
        if (!contourParameter || !contourParameter.values) return -1;
        if (Array.isArray(category)) {
            return category.map(c => contourParameter.values.indexOf(c));
        }
        return contourParameter.values.indexOf(category);
    }, [contourParameter, category]);
    const filteredIndices = useMemo(() => {
        if (categoryValueIndex === -1) return [];
        if (Array.isArray(categoryValueIndex)) {
            return data.filter(i => categoryValueIndex.includes(contourParameter.data[i]));
        }
        return data.filter(i => contourParameter.data[i] === categoryValueIndex);
    }, [data, categoryValueIndex, contourParameter]);
    return filteredIndices;
}



/**
 * This assumes that the current chart context has a `config.region` key that refers to a region with `viv_image` in the data store.
 */
export function useImgUrl(): string {
    const region = useRegion();
    const avivator = useDataStore().regions.avivator;
    const url = useMemo(() => {
        // if (config.imageURL) return config.imageURL; //deprecated
        const i = region.viv_image;
        if (avivator.base_url?.startsWith("http")) {
            const url = new URL(i.url || i.file, avivator.base_url);
            return url.href;
        }
        // todo - this is a bit of a mess, and should be cleaned up / better documented.
        const url = i.url ? i.url : getProjectURL(avivator.base_url) + i.file;
        return url;
    }, [region, avivator]);
    return url;
}

/**
 * Returns a list of `dataSources` that are currently available.
 * This is `observable` and will update when new data sources are added etc.
 */
export function useDataSources() {
    return window.mdv.chartManager?.dataSources;
}
