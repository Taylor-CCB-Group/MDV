import { readZarr, type SpatialData } from "@spatialdata/core";

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

export type SpatialDataReader = {
    /** The read of `url`, starting one if there is none. Safe to call during render. */
    read(url: string): SpatialDataRead;
    /**
     * Starts any failed read again. For a user action only: retrying from render
     * would loop with the render that asked for it.
     */
    retry(): void;
    /** Called whenever a read settles or is replaced by `retry`. */
    subscribe(listener: () => void): () => void;
};

/**
 * One `readZarr` per store URL, held by whoever owns the reader.
 *
 * Opening a store is not cheap. On a v2 store zarrextra asks for `zarr.json` before
 * `.zarray` for every array, and zarrita asks for `.zattrs` on every node. The
 * consolidated store sends whatever it does not know to the network and does not
 * remember misses, so each of those is a sequential 404. Each upstream
 * `SpatialDataProvider` does its own read, so the layer dialog — a separate React
 * tree — used to repeat all of it on every open (#610).
 *
 * The chart owns one of these and both its trees read through it. Reads are kept
 * per URL, not in a single slot: while the trees disagree about the URL (one has
 * re-rendered after a region change and the other has not) each would otherwise
 * evict the other's read on every notification. A chart sees very few URLs, so
 * nothing is evicted. A failed read stays failed until `retry`.
 */
export function createSpatialDataReader(load: (url: string) => Promise<SpatialData> = readZarr): SpatialDataReader {
    const reads = new Map<string, SpatialDataRead>();
    const listeners = new Set<() => void>();
    const notify = () => {
        for (const listener of listeners) listener();
    };
    const start = (url: string): SpatialDataRead => {
        const entry: SpatialDataRead = {
            url,
            promise: load(url),
            state: { spatialData: null, loading: true, error: null },
        };
        entry.promise.then(
            (spatialData) => {
                entry.state = { spatialData, loading: false, error: null };
                notify();
            },
            (error) => {
                entry.state = {
                    spatialData: null,
                    loading: false,
                    error: error instanceof Error ? error : new Error(String(error)),
                };
                notify();
            },
        );
        reads.set(url, entry);
        return entry;
    };
    return {
        read: (url) => reads.get(url) ?? start(url),
        retry: () => {
            const failed = [...reads.values()].filter((entry) => entry.state.error);
            if (!failed.length) return;
            for (const entry of failed) start(entry.url);
            notify();
        },
        subscribe: (listener) => {
            listeners.add(listener);
            return () => {
                listeners.delete(listener);
            };
        },
    };
}
