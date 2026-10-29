import { type SpatialData, readZarr } from "@spatialdata/core";

/** Same shape as `useSpatialData` from `@spatialdata/react`, which this replaces in MDV. */
export type SpatialDataState = {
    spatialData: SpatialData | null;
    loading: boolean;
    error: Error | null;
};

export type SpatialDataRead = {
    readonly url: string;
    readonly promise: Promise<SpatialData>;
    /** Replaced (not mutated) when the read settles, so its identity is a change signal. */
    state: SpatialDataState;
};

export const NO_SPATIAL_DATA: SpatialDataState = { spatialData: null, loading: false, error: null };

/**
 * One `readZarr` per store URL, held by whoever owns the returned function.
 *
 * Opening a store is not cheap: zarrextra walks every node and, on a v2 store, asks
 * for `zarr.json` before `.zarray`/`.zgroup`. That key is never in v2 consolidated
 * metadata and zarrita does not remember misses, so every node is a sequential
 * network 404. Each upstream `SpatialDataProvider` does its own read, so the layer
 * dialog — a separate React tree — used to repeat all of it on every open (#610).
 *
 * The chart owns one of these and both its trees read through it. A failed read
 * stays failed until the URL changes; retrying here would loop with the render
 * that asked for it.
 */
export function createSpatialDataReader(read: (url: string) => Promise<SpatialData> = readZarr) {
    let current: SpatialDataRead | undefined;
    return (url: string): SpatialDataRead => {
        if (current?.url === url) return current;
        const entry: SpatialDataRead = {
            url,
            promise: read(url),
            state: { spatialData: null, loading: true, error: null },
        };
        entry.promise.then(
            (spatialData) => {
                entry.state = { spatialData, loading: false, error: null };
            },
            (error) => {
                entry.state = {
                    spatialData: null,
                    loading: false,
                    error: error instanceof Error ? error : new Error(String(error)),
                };
            },
        );
        current = entry;
        return entry;
    };
}
