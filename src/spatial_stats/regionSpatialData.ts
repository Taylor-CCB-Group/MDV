import type { SpatialData } from "@spatialdata/core";
import type DataStore from "@/datastore/DataStore";
import { spatialDataUrlForRegion } from "@/react/spatialdata/shared_spatial_data";
import { createSpatialDataReader } from "@/react/spatialdata/spatial_data_reader";
import type { ConvertedRegion } from "./spatialDataFrame";

/** Store reads for spatial statistics outside any chart (a chart keeps its own reader). */
const reader = createSpatialDataReader();

/** The SpatialData regions of a datastore: entries with a store behind them. Converted projects also
 *  carry image-only regions (e.g. a lossy OME-TIFF), which don't count. */
export function spatialRegionKeys(dataStore: DataStore): string[] {
    const all: Record<string, ConvertedRegion> = dataStore.regions?.all_regions ?? {};
    return Object.keys(all).filter((k) => spatialDataUrlForRegion(all[k]));
}

export type RegionStore = {
    key: string | undefined;
    entry: ConvertedRegion | undefined;
    /** Absent when the region has no store, or reading it failed (`error`). */
    sdata?: SpatialData;
    error?: string;
};

/** `region`, or the only SpatialData region when none is given, and its store (read once per URL). */
export async function readRegionStore(dataStore: DataStore, region?: string): Promise<RegionStore> {
    const keys = spatialRegionKeys(dataStore);
    const key = region ?? (keys.length === 1 ? keys[0] : undefined);
    const entry: ConvertedRegion | undefined = key ? dataStore.regions?.all_regions[key] : undefined;
    const url = entry ? spatialDataUrlForRegion(entry) : null;
    if (!url) return { key, entry };
    try {
        return { key, entry, sdata: await reader.read(url).promise };
    } catch (e) {
        return { key, entry, error: `store read failed: ${e}` };
    }
}
