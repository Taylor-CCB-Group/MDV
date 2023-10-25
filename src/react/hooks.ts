import { useEffect, useMemo, useState } from "react";
import BaseChart from "../charts/BaseChart";
import { loadOmeTiff, getChannelStats } from "@hms-dbmi/viv";
import { BaseReactChart } from "./components/BaseReactChart";

export function useChartSize(chart: BaseChart) {
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
 * consider context API for some of these things?
 */
export function useChartID(chart: BaseChart) {
    return chart.config.id;
}

export function useFilteredIndices(chart: BaseReactChart<any>) {
    const filterArray = chart.dataStore.filterArray as Uint8Array;
    const [filteredIndices, setFilteredIndices] = useState<Uint32Array>(new Uint32Array(filterArray.length).map((_, i) => i));
    useEffect(() => {
        chart.dataStore.getFilteredIndices().then(setFilteredIndices);
    }, [chart.dataStore._filteredIndicesPromise]);
    return filteredIndices;
}

export function useParamColumns(chart: BaseChart) {
    const { columnIndex } = chart.dataStore;
    const columns = useMemo(() => chart.config.param.map(name => columnIndex[name]), [chart.config.param]);
    return columns;
}

type OME_TIFF = Awaited<ReturnType<typeof loadOmeTiff>>;
export function useOmeTiff(url: string) {
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

// mobx for everything?
// we could have a config.channelsState, in which case it shouldn't need
// too much more architecture to make it work.
// Color
export type ChannelsState = {
    channelsVisible: boolean[],
    contrastLimits: [number, number][],
    colors: [number, number, number][],
    domains: [number, number][],
    selections: { z: number, c: number, t: number }[],
    ids: string[]
}
const DEFAUlT_CHANNEL_STATE = {
    channelsVisible: [],
    contrastLimits: [],
    colors: [],
    domains: [],
    selections: [],
    ids: [],
    loader: [{ labels: [], shape: [] }],
    image: 0
};
