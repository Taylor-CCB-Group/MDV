import { useEffect, useLayoutEffect, useMemo, useState } from "react";
import { loadOmeTiff, getChannelStats } from "@hms-dbmi/viv";
import { useChart, useDataStore } from "./context";
import type { OME_TIFF } from "./components/avivatorish/state";
import { getProjectURL } from "../dataloaders/DataLoaderUtil";
import { getRandomString } from "../utilities/Utilities";
import { action } from "mobx";
import type { DataColumn } from "../charts/charts";

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
        return chart.config.param.map(name => columnIndex[name])
    }, [chart.config.param, columnIndex]);
    return columns;
}

// slightly rough - we might have multiple images in a config, or generally think about this differently
// this is breaking rules of hooks etc in short term while I figure out how to do this properly
export function useImgUrl(): string {
    const config = useConfig() as any;
    if (config.imageURL) return config.imageURL;
    // see VivScatterPlot.afterAppCreation() ...
    const { regions } = useDataStore();
    if (!regions) {
        //throw `No image URL provided and no regions found in data store`; // plenty of other ways this could go wrong}
        console.warn("No image URL provided and no regions found in data store"); // plenty of other ways this could go wrong}
        return '';
    }
    const i = regions.all_regions[config.region].viv_image;
    return i.url ? i.url : getProjectURL(regions.avivator.base_url) + i.file;
}
