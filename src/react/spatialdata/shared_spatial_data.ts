import { useCallback, useSyncExternalStore } from "react";
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
 * first render, and every tree sees a `retry`.
 */
export function useChartSpatialData(): SpatialDataState {
    const reader = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>().spatialDataReader;
    const url = spatialDataUrlForRegion(useRegion());
    const getSnapshot = useCallback(() => (url ? reader.read(url).state : NO_SPATIAL_DATA), [reader, url]);
    return useSyncExternalStore(reader.subscribe, getSnapshot);
}
