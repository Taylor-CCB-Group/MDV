import { useEffect, useReducer } from "react";
import { getProjectURL } from "@/dataloaders/DataLoaderUtil";
import { useChart } from "../context";
import { useRegion } from "../hooks";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "../components/SpatialDataMDVReact";
import { NO_SPATIAL_DATA, type SpatialDataState } from "./spatial_data_reader";

export function spatialDataUrlForRegion(region: unknown): string | null {
    if (!region || typeof region !== "object" || !("spatial" in region)) return null;
    const spatial = region.spatial;
    if (!spatial || typeof spatial !== "object" || !("file" in spatial)) return null;
    const file = spatial.file;
    if (typeof file !== "string" || !file) return null;
    return getProjectURL(`spatial/${file}`);
}

/**
 * The chart's `SpatialData`, for any tree under its `ChartProvider` — the canvas and
 * the layer dialog alike. A tree that mounts after the read settled gets it on its
 * first render.
 */
export function useChartSpatialData(): SpatialDataState {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const url = spatialDataUrlForRegion(useRegion());
    const read = url ? chart.readSpatialData(url) : null;
    const [, rerender] = useReducer((n: number) => n + 1, 0);
    const pending = read?.state.loading ?? false;
    useEffect(() => {
        if (!read || !pending) return;
        let mounted = true;
        // the reader's own handlers were attached first, so `read.state` is final here
        const settle = () => mounted && rerender();
        read.promise.then(settle, settle);
        return () => {
            mounted = false;
        };
    }, [read, pending]);
    return read?.state ?? NO_SPATIAL_DATA;
}
