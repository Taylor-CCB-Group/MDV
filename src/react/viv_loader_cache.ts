import { createLoader, type UrlOrFiles } from "./components/avivatorish/utils";

type VivLoaderResult = Awaited<ReturnType<typeof createLoader>>;

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
    for (const [key] of removable) {
        if (toRemove <= 0) break;
        const current = vivLoaderCache.get(key);
        if (!current || current.promise) continue;
        vivLoaderCache.delete(key);
        toRemove -= 1;
    }
}

export async function getOrCreateVivLoader(
    urlOrFile: UrlOrFiles,
    handleOffsetsNotFound: (arg0: boolean) => void,
    handleLoaderError: (msg: string | null) => void,
) {
    const cacheKey = getCacheKey(urlOrFile);
    if (!cacheKey) {
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
            return existing.value;
        }
        if (existing.promise) {
            existing.subscribers.add(subscriber);
            try {
                return await existing.promise;
            } finally {
                existing.subscribers.delete(subscriber);
            }
        }
    }

    const entry: CacheEntry = {
        lastAccessMs: Date.now(),
        subscribers: new Set([subscriber]),
    };
    const promise = createLoader(
        urlOrFile,
        (value) => notifyOffsets(entry, value),
        (message) => notifyLoaderError(entry, message),
    )
        .then((result) => {
            entry.value = result;
            entry.promise = undefined;
            entry.subscribers.clear();
            touchEntry(entry);
            pruneVivLoaderCache();
            return result;
        })
        .catch((error) => {
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

