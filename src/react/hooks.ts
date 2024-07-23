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

/**
 * **Should only be used in OmeTiffProvider.**
 * 
 * Get an OME tiff from the chart's config. As of now, this'll either look for config.imageURL,
 * or attempt to ascertain a URL based on regions & other config...
 * 
 * This is likely to change in future - we want to support other image formats, and likely other forms of config etc.
 */
export function useOmeTiffLoader() {
    //uh oh... different ome state for every use of this hook... maybe we should use a context?
    //or a zustand store?
    //if we're using multiple images in the same Deck thing, may not fit well with context.
    //(also wouldn't make sense for no-args function, so we could for now...)
    ///--- for now, I'm going to reduce the places I call useOmeTiff()
    const url = useImgUrl();
    const [tiff, setTiff] = useState<OME_TIFF>();
    useEffect(() => {
        if (url.endsWith('.ome.tif') || url.endsWith('.ome.tiff')) {
            loadOmeTiff(url).then(setTiff);
        }
        else {
            throw `Only .ome.tif currently supported - ${url}`;
            // const _url = new URL(url, document.baseURI).href;
            // loadOmeZarr(_url, { type: 'multiscales' }).then(setTiff);
        }
    }, [url]);
    return tiff;
}

/**
 * @deprecated in favour of avivatorish hooks for now.
 * 
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

