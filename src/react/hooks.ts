import { useEffect, useMemo, useState } from "react";
import BaseChart from "../charts/BaseChart";
import { DataModel } from "../table/DataModel";
import { loadOmeTiff, getChannelStats } from "@hms-dbmi/viv";

export function useChartSize(chart: BaseChart) {
    const div = chart.contentDiv;
    const [width, setWidth] = useState(div.clientWidth);
    const [height, setHeight] = useState(div.clientHeight);
    useEffect(() => {
        const resize = () => {
            setWidth(div.clientWidth);
            setHeight(div.clientHeight);
        };
        const observer = new ResizeObserver(resize);
        observer.observe(div);
        return () => observer.unobserve(div);
    }, [div]);
    return [width, height];
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

/**
 * This hook allows you to use a named member from a chart's config object as a react state variable.
 * 
 * ## notes:
 * *This design is not final and may actually be more broken than it appears.*
 * Probably better to use mobx instead, at least in the short term.
 * 
 * @param chart - the chart to get the config item from. If the class of the chart explicity declares the config type,
 * this will be used to infer the type of the returned value.
 * 
 * Otherwise, it will be `any`.
 * @param key 
 */
export function useConfigItem<T extends BaseChart, K extends keyof T["config"]>(chart: T, key: K) {
    const [value, setValue] = useState<T["config"][K]>(chart.config[key]);
    useEffect(() => {
        Object.defineProperty(chart.config, key, {
            get() {
                // **this is closed on the initial value when the property is defined; not what we want.**
                return value;
            },
            set(val) {
                setValue(val);
            }
        });
    }, [chart.config, key, chart.config[key]]);
    return value;
}