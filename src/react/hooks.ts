import { useEffect, useMemo, useState } from "react";
import BaseChart from "../charts/BaseChart";
import { DataModel } from "../table/DataModel";
import { loadOmeTiff, getChannelStats } from "@hms-dbmi/viv";

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

export function useDataModel(chart: BaseChart) {
    const { dataStore } = chart;
    const [dataModel, setDataModel] = useState(() => new DataModel(dataStore, { autoUpdate: true }));
    useEffect(() => {
        dataModel.setColumns(chart.config.param);
        dataModel.updateModel();
    }, [dataModel]);
    // as of now, config.param can't change, but I think it should be able to.
    // this will need to be tested in that case.
    useEffect(() => {
        dataModel.setColumns(chart.config.param);
        dataModel.updateModel();
    }, [chart.config.param]);
    return dataModel;
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
