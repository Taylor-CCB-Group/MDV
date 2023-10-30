import { useEffect, useMemo, useState } from "react";
import { loadOmeTiff, getChannelStats } from "@hms-dbmi/viv";
import { BaseConfig } from "./components/BaseReactChart";
import { useChart, useDataStore } from "./context";
import { ChannelsState, VivConfig } from "./viv_state";
import { VivRoiConfig } from "./components/VivMDVReact";
import { getProjectURL } from "../dataloaders/DataLoaderUtil";

/**
 * Get the chart's config.
 * 
 * Must be used within a ChartContext.
 * 
 * Provided type parameter is not checked - in future it could probably
 * be inferred from the chart type.
 */
export function useConfig<T extends BaseConfig>() {
    const { config } = useChart();
    return config as T; //todo: strict/inferred typing
}

export function useChartSize() {
    const chart = useChart();
    // return chart.config.size; // not so well behaved?
    const div = chart.contentDiv;
    const [size, setSize] = useState([div.clientWidth, div.clientHeight]);
    useEffect(() => {
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
export function useChartID() {
    const chart = useChart();
    return chart.config.id;
}

/**
 * Get a {Uint32Array} of the currently filtered indices.
 * When the selection changes, this will asynchronously update.
 * All users of the same data store share a reference to the same array.
 */
export function useFilteredIndices() {
    const dataStore = useDataStore();
    const [filteredIndices, setFilteredIndices] = useState(new Uint32Array());
    useEffect(() => {
        dataStore.getFilteredIndices().then(setFilteredIndices);
        // using _filteredIndicesPromise as a dependency is working reasonably well,
        // but possibly needs a bit more thought.
    }, [dataStore._filteredIndicesPromise]);
    return filteredIndices;
}

export function useParamColumns() {
    const chart = useChart();
    const { columnIndex } = chart.dataStore;
    const columns = useMemo(() => chart.config.param.map(name => columnIndex[name]), [chart.config.param]);
    return columns;
}

type OME_TIFF = Awaited<ReturnType<typeof loadOmeTiff>>;

// slightly rough - we might have multiple images in a config, or generally think about this differently
// this is breaking rules of hooks etc in short term while I figure out how to do this properly
function useImgUrl() {
    const config = useConfig() as any;
    if (config.imageURL) return config.imageURL;
    // see VivScatterPlot.afterAppCreation() ...
    const { regions } = useDataStore();
    if (!regions) throw `No image URL provided and no regions found in data store`; // plenty of other ways this could go wrong
    const i = regions.all_regions[config.region].viv_image;
    return i.url ? i.url : getProjectURL(regions.avivator.base_url) + i.file;
}

/**
 * Get an OME tiff from the chart's config. As of now, this'll either look for config.imageURL,
 * or attempt to ascertain a URL based on regions & other config...
 * 
 * This is likely to change in future - we want to support other image formats, and likely other forms of config etc.
 */
export function useOmeTiff() {
    const url = useImgUrl();
    const [tiff, setTiff] = useState<OME_TIFF>();
    useEffect(() => {
        loadOmeTiff(url).then(setTiff);
    }, [url]);
    return tiff;
}

/**
 * Get channel statistics for a given channel of an OME tiff.
 * @param ome 
 * @param channel - channel index. In future we should support selections with z/c/t.
 */
export function useChannelStats(ome: OME_TIFF, channel: number) {
    const [channelStats, setChannelStats] = useState<ReturnType<typeof getChannelStats> | undefined>();
    useEffect(() => {
        if (!ome) return;
        const loader = ome.data[ome.data.length - 1];
        const selection = { z: 0, c: channel, t: 0 };
        loader.getRaster({ selection }).then((raster) => {
            setChannelStats(getChannelStats(raster.data));
        });
    }, [ome, channel]);
    return channelStats;
}


// -- signature/???
export function useChannelsState(state: ChannelsState) {
    const [channelsState, setChannelsState] = useState(state);
    return channelsState;
}