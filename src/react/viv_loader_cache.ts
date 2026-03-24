import { observable, runInAction } from "mobx";
import { createLoader, type UrlOrFiles } from "./components/avivatorish/utils";

const LOG_PREFIX = "[vivLoaderCache]";
/** Log tile loads that exceed this duration (ms). Set to Infinity to disable. */
const SLOW_TILE_LOAD_LOG_MS = 750;

type VivLoaderResult = Awaited<ReturnType<typeof createLoader>>;
type PixelData = { data: ArrayLike<number>; width: number; height: number };
type TileArgs = { x: number; y: number; selection: unknown; signal?: AbortSignal };
type RasterArgs = { selection: unknown; signal?: AbortSignal };
type PixelSourceLike = {
    getTile?: (args: TileArgs) => Promise<PixelData>;
    getRaster?: (args: RasterArgs) => Promise<PixelData>;
    pool?: unknown;
};

type LoaderCallbacks = {
    handleOffsetsNotFound: (arg0: boolean) => void;
    handleLoaderError: (msg: string | null) => void;
};

type CacheEntry = {
    lastAccessMs: number;
    value?: VivLoaderResult;
    promise?: Promise<VivLoaderResult>;
    subscribers: Set<LoaderCallbacks>;
};

const DEFAULT_MAX_CACHE_ENTRIES = 12;
const vivLoaderCache = new Map<string, CacheEntry>();
const DEFAULT_MAX_TILE_CACHE_ENTRIES = 800;

type TileCacheEntry = {
    lastAccessMs: number;
    value?: PixelData;
    promise?: Promise<PixelData>;
};
const vivTileCache = new Map<string, TileCacheEntry>();

const telemetry = {
    /** Resolved loader returned from cache (warm entry). */
    loaderCacheHitValue: 0,
    /** Awaited in-flight loader promise for the same URL. */
    loaderCacheHitInflight: 0,
    /** Started a new createLoader for this cache key. */
    loaderCacheMissNew: 0,
    /** URL/file could not be keyed; bypassed cache. */
    loaderPassthroughNonStringUrl: 0,
    /** createLoader promise rejected after cache insert. */
    loaderCreateFailures: 0,
    /** Entries removed by pruneVivLoaderCache. */
    loaderEntriesPruned: 0,

    tileCacheHit: 0,
    tileCacheInflightReuse: 0,
    tileCacheMissLoad: 0,
    tileLoadFailures: 0,
    tileEntriesPruned: 0,

    pixelSourcesWrapped: 0,

    /** Sum of createLoader wall time for completed loads (ms). */
    loaderCreateDurationMsTotal: 0,
    loaderCreateCompletedCount: 0,
};

function truncateForLog(s: string, max = 160) {
    if (s.length <= max) return s;
    return `${s.slice(0, max)}…(${s.length} chars)`;
}

export type VivLoaderCacheTelemetry = Readonly<typeof telemetry> & {
    vivLoaderCacheSize: number;
    vivTileCacheSize: number;
    /** Mean createLoader duration for completed loads (ms); 0 if none. */
    avgLoaderCreateMs: number;
    /**
     * Fraction of tile/raster requests served from warm cache or shared in-flight work
     * (vs cold decode). 0 if no tile requests recorded.
     */
    tileLayerReuseRatio: number;
};

/** Counters and timings for the Viv loader / tile caches (devtools, tests, support). */
export function getVivLoaderCacheTelemetry(): VivLoaderCacheTelemetry {
    const avgLoaderCreateMs =
        telemetry.loaderCreateCompletedCount > 0
            ? telemetry.loaderCreateDurationMsTotal / telemetry.loaderCreateCompletedCount
            : 0;
    const tileTotal =
        telemetry.tileCacheHit + telemetry.tileCacheInflightReuse + telemetry.tileCacheMissLoad;
    const tileLayerReuseRatio =
        tileTotal > 0 ? (telemetry.tileCacheHit + telemetry.tileCacheInflightReuse) / tileTotal : 0;
    return {
        ...telemetry,
        vivLoaderCacheSize: vivLoaderCache.size,
        vivTileCacheSize: vivTileCache.size,
        avgLoaderCreateMs,
        tileLayerReuseRatio,
    };
}

/**
 * Live snapshot for debug UIs (e.g. Debug JSON dialog). Updated on a short throttle so tile-heavy
 * panning does not force excessive re-renders.
 */
export const vivLoaderCacheTelemetryObservable = observable({
    snapshot: getVivLoaderCacheTelemetry(),
});

const TELEMETRY_PUBLISH_MS = 250;
let telemetryPublishTimer: ReturnType<typeof setTimeout> | undefined;

function publishVivLoaderTelemetry() {
    if (telemetryPublishTimer !== undefined) return;
    telemetryPublishTimer = setTimeout(() => {
        telemetryPublishTimer = undefined;
        runInAction(() => {
            vivLoaderCacheTelemetryObservable.snapshot = getVivLoaderCacheTelemetry();
        });
    }, TELEMETRY_PUBLISH_MS);
}

const pixelSourceWrapperCache = new WeakMap<object, object>();

function touchEntry(entry: CacheEntry) {
    entry.lastAccessMs = Date.now();
}

function notifyOffsets(entry: CacheEntry, value: boolean) {
    for (const subscriber of entry.subscribers) {
        subscriber.handleOffsetsNotFound(value);
    }
}

function notifyLoaderError(entry: CacheEntry, message: string | null) {
    for (const subscriber of entry.subscribers) {
        subscriber.handleLoaderError(message);
    }
}

function getCacheKey(urlOrFile: UrlOrFiles) {
    // Conservative first pass: cache only string URLs.
    return typeof urlOrFile === "string" ? urlOrFile : null;
}

export function pruneVivLoaderCache(maxEntries = DEFAULT_MAX_CACHE_ENTRIES) {
    if (vivLoaderCache.size <= maxEntries) return;

    const removable = [...vivLoaderCache.entries()]
        .filter(([, entry]) => !entry.promise)
        .sort((a, b) => a[1].lastAccessMs - b[1].lastAccessMs);

    let toRemove = vivLoaderCache.size - maxEntries;
    let removed = 0;
    for (const [key] of removable) {
        if (toRemove <= 0) break;
        const current = vivLoaderCache.get(key);
        if (!current || current.promise) continue;
        vivLoaderCache.delete(key);
        toRemove -= 1;
        removed += 1;
    }
    if (removed > 0) {
        telemetry.loaderEntriesPruned += removed;
        console.debug(
            `${LOG_PREFIX} pruned ${removed} loader entries (max=${maxEntries}, remaining=${vivLoaderCache.size})`,
        );
        publishVivLoaderTelemetry();
    }
}

function stableSerialize(value: unknown): string {
    if (value === null) return "null";
    if (value === undefined) return "undefined";
    const t = typeof value;
    if (t === "number" || t === "boolean" || t === "string") return JSON.stringify(value);
    if (Array.isArray(value)) return `[${value.map((v) => stableSerialize(v)).join(",")}]`;
    if (t === "object") {
        const obj = value as Record<string, unknown>;
        const keys = Object.keys(obj).sort();
        return `{${keys.map((k) => `${JSON.stringify(k)}:${stableSerialize(obj[k])}`).join(",")}}`;
    }
    return JSON.stringify(String(value));
}

function pruneVivTileCache(maxEntries = DEFAULT_MAX_TILE_CACHE_ENTRIES) {
    if (vivTileCache.size <= maxEntries) return;
    const removable = [...vivTileCache.entries()]
        .filter(([, entry]) => !entry.promise)
        .sort((a, b) => a[1].lastAccessMs - b[1].lastAccessMs);
    let toRemove = vivTileCache.size - maxEntries;
    let removed = 0;
    for (const [key] of removable) {
        if (toRemove <= 0) break;
        const current = vivTileCache.get(key);
        if (!current || current.promise) continue;
        vivTileCache.delete(key);
        toRemove -= 1;
        removed += 1;
    }
    if (removed > 0) {
        telemetry.tileEntriesPruned += removed;
        console.debug(
            `${LOG_PREFIX} pruned ${removed} tile entries (max=${maxEntries}, remaining=${vivTileCache.size})`,
        );
        publishVivLoaderTelemetry();
    }
}

async function getOrCreateTileData(
    key: string,
    signal: AbortSignal | undefined,
    load: () => Promise<PixelData>,
) {
    const existing = vivTileCache.get(key);
    if (existing?.value) {
        telemetry.tileCacheHit += 1;
        publishVivLoaderTelemetry();
        existing.lastAccessMs = Date.now();
        return existing.value;
    }
    // If caller provides a signal, avoid sharing in-flight work keyed to that signal.
    if (!signal && existing?.promise) {
        telemetry.tileCacheInflightReuse += 1;
        publishVivLoaderTelemetry();
        existing.lastAccessMs = Date.now();
        return existing.promise;
    }

    const entry: TileCacheEntry = existing ?? { lastAccessMs: Date.now() };
    entry.lastAccessMs = Date.now();
    const timedLoad = () => {
        telemetry.tileCacheMissLoad += 1;
        publishVivLoaderTelemetry();
        const t0 = performance.now();
        return load().then((value) => {
            const dt = performance.now() - t0;
            if (dt >= SLOW_TILE_LOAD_LOG_MS) {
                console.warn(
                    `${LOG_PREFIX} slow tile load ${dt.toFixed(0)}ms`,
                    truncateForLog(key),
                );
            }
            return value;
        });
    };

    const promise = timedLoad()
        .then((value) => {
            entry.value = value;
            entry.promise = undefined;
            entry.lastAccessMs = Date.now();
            vivTileCache.set(key, entry);
            pruneVivTileCache();
            publishVivLoaderTelemetry();
            return value;
        })
        .catch((error) => {
            telemetry.tileLoadFailures += 1;
            publishVivLoaderTelemetry();
            console.warn(`${LOG_PREFIX} tile load failed`, truncateForLog(key), error);
            const current = vivTileCache.get(key);
            if (current?.promise === promise) {
                vivTileCache.delete(key);
            }
            throw error;
        });

    if (!signal) {
        entry.promise = promise;
        vivTileCache.set(key, entry);
        pruneVivTileCache();
    }

    return promise;
}

function isLikelyZarrPixelSource(source: PixelSourceLike) {
    // TiffPixelSource exposes a `pool` property; keep wrapping conservative to avoid side effects.
    return typeof source.getTile === "function" && typeof source.getRaster === "function" && !("pool" in source);
}

function wrapPixelSource(source: PixelSourceLike, sourceKey: string, levelKey: string) {
    if (!isLikelyZarrPixelSource(source)) return source;
    const existing = pixelSourceWrapperCache.get(source as object);
    if (existing) return existing as PixelSourceLike;

    const wrapped = new Proxy(source as object, {
        get(target, prop, receiver) {
            if (prop === "getTile") {
                return async (args: TileArgs) => {
                    const key = `${sourceKey}|${levelKey}|tile|x:${args.x}|y:${args.y}|sel:${stableSerialize(args.selection)}`;
                    return getOrCreateTileData(
                        key,
                        args.signal,
                        () => (source.getTile as NonNullable<PixelSourceLike["getTile"]>)(args),
                    );
                };
            }
            if (prop === "getRaster") {
                return async (args: RasterArgs) => {
                    const key = `${sourceKey}|${levelKey}|raster|sel:${stableSerialize(args.selection)}`;
                    return getOrCreateTileData(
                        key,
                        args.signal,
                        () => (source.getRaster as NonNullable<PixelSourceLike["getRaster"]>)(args),
                    );
                };
            }
            const value = Reflect.get(target, prop, receiver);
            return typeof value === "function" ? value.bind(target) : value;
        },
    });
    pixelSourceWrapperCache.set(source as object, wrapped);
    telemetry.pixelSourcesWrapped += 1;
    publishVivLoaderTelemetry();
    console.debug(`${LOG_PREFIX} wrapped Zarr-like PixelSource for tile cache`, truncateForLog(sourceKey));
    return wrapped as PixelSourceLike;
}

function wrapLoaderResultForTileCache(result: VivLoaderResult, sourceKey: string): VivLoaderResult {
    if (Array.isArray(result)) {
        return result.map((image, imageIndex) => ({
            ...image,
            data: image.data.map((source, levelIndex) =>
                wrapPixelSource(source as PixelSourceLike, sourceKey, `img:${imageIndex}|lvl:${levelIndex}`),
            ),
        })) as VivLoaderResult;
    }
    if (result && typeof result === "object" && "data" in result && Array.isArray(result.data)) {
        return {
            ...result,
            data: result.data.map((source, levelIndex) =>
                wrapPixelSource(source as PixelSourceLike, sourceKey, `img:0|lvl:${levelIndex}`),
            ),
        } as VivLoaderResult;
    }
    return result;
}

export async function getOrCreateVivLoader(
    urlOrFile: UrlOrFiles,
    handleOffsetsNotFound: (arg0: boolean) => void,
    handleLoaderError: (msg: string | null) => void,
) {
    const cacheKey = getCacheKey(urlOrFile);
    if (!cacheKey) {
        telemetry.loaderPassthroughNonStringUrl += 1;
        publishVivLoaderTelemetry();
        console.debug(`${LOG_PREFIX} uncached createLoader (non-string UrlOrFiles)`);
        return createLoader(urlOrFile, handleOffsetsNotFound, handleLoaderError);
    }

    const subscriber: LoaderCallbacks = {
        handleOffsetsNotFound,
        handleLoaderError,
    };

    const existing = vivLoaderCache.get(cacheKey);
    if (existing) {
        touchEntry(existing);
        if (existing.value) {
            telemetry.loaderCacheHitValue += 1;
            publishVivLoaderTelemetry();
            console.debug(`${LOG_PREFIX} loader cache hit (ready)`, truncateForLog(cacheKey));
            return existing.value;
        }
        if (existing.promise) {
            telemetry.loaderCacheHitInflight += 1;
            publishVivLoaderTelemetry();
            console.debug(`${LOG_PREFIX} loader cache hit (in-flight)`, truncateForLog(cacheKey));
            existing.subscribers.add(subscriber);
            try {
                return await existing.promise;
            } finally {
                existing.subscribers.delete(subscriber);
            }
        }
    }

    telemetry.loaderCacheMissNew += 1;
    publishVivLoaderTelemetry();
    console.debug(`${LOG_PREFIX} loader cache miss, creating`, truncateForLog(cacheKey));

    const entry: CacheEntry = {
        lastAccessMs: Date.now(),
        subscribers: new Set([subscriber]),
    };
    const loadStartedAt = performance.now();
    const promise = createLoader(
        urlOrFile,
        (value) => notifyOffsets(entry, value),
        (message) => notifyLoaderError(entry, message),
    )
        .then((result) => {
            const dt = performance.now() - loadStartedAt;
            telemetry.loaderCreateDurationMsTotal += dt;
            telemetry.loaderCreateCompletedCount += 1;
            publishVivLoaderTelemetry();
            console.debug(
                `${LOG_PREFIX} createLoader finished in ${dt.toFixed(0)}ms`,
                truncateForLog(cacheKey),
            );
            const wrapped = wrapLoaderResultForTileCache(result, cacheKey);
            entry.value = wrapped;
            entry.promise = undefined;
            entry.subscribers.clear();
            touchEntry(entry);
            pruneVivLoaderCache();
            return wrapped;
        })
        .catch((error) => {
            telemetry.loaderCreateFailures += 1;
            publishVivLoaderTelemetry();
            console.warn(`${LOG_PREFIX} createLoader failed`, truncateForLog(cacheKey), error);
            vivLoaderCache.delete(cacheKey);
            entry.subscribers.clear();
            throw error;
        });

    entry.promise = promise;
    vivLoaderCache.set(cacheKey, entry);
    pruneVivLoaderCache();

    try {
        return await promise;
    } finally {
        entry.subscribers.delete(subscriber);
    }
}

