import { useEffect, useMemo, useState } from "react";
import BaseChart from "../charts/BaseChart";
import { DataModel } from "../table/DataModel";
import { loadOmeTiff } from "@hms-dbmi/viv";

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

export function useOmeTiff(url: string) {
    const [tiff, setTiff] = useState<Awaited<ReturnType<typeof loadOmeTiff>>>();
    useEffect(() => {
        loadOmeTiff(url).then(setTiff);
    }, [url]);
    return tiff;
}

