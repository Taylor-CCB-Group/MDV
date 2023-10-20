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

type OME = Awaited<ReturnType<typeof loadOmeTiff>>;
export function useOmeTiff(url: string) {
    const [tiff, setTiff] = useState<OME>();
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
export function useChannelStats(ome: OME, channel: number) {
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
 * Bridge between MDV settings and react state... there is a lot more to think about longer-term.
 * I don't much like this approach, but it's a start.
 * Could consider MobX or something similar.
 * Also consider a 'useSettings' hook which will call `chart.getSettings()` and process the result.
 * @param chart
 * @param key 
 */
export function useConfigItem(chart: BaseChart, key: string) {
    const [value, setValue] = useState(chart.config[key]);
    useEffect(() => {
        // test cases to consider:
        // - referring to the same chart.config[key] in multiple places
        //   how bad will it be repeatedly calling defineProperty?
        // - Under what other circumstances may we be liable to call defineProperty multiple times?
        Object.defineProperty(chart.config, key, {
            get() {
                return value;
            },
            set(val) {
                setValue(val);
            }
        });
    }, [chart.config, key, chart.config[key]]);
    return value;
}