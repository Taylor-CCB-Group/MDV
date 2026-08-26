import { ensureWorkers } from "@spatialdata/vis";

// Vite bundles the core parquet worker and hands back a runtime URL. The import is
// what does the work — it asks Vite to *build* the worker, shared chunks and its own
// parquet-wasm included. The worker cannot be shipped self-contained, so there is no
// default URL that survives being re-bundled into MDV's `assets/`; see
// "Bundling into an application" in the SpatialData.js docs.
import parquetWorkerUrl from "@spatialdata/core/parquet-worker?worker&url";

/**
 * Enable the core parquet worker before any spatial parquet work.
 *
 * Not optional for anything past the first capped read. Core's `defaultEnabled` is
 * `false`, and `loadPointsMatchingFeatureCodes` — the feature-index scan that fetches a
 * selected feature's points from beyond the resident window — throws outright without
 * it. Every other caller degrades to a main-thread decode instead, so a worker that
 * fails to load costs performance rather than correctness.
 *
 * It also moves the heavy decodes off the main thread: the codes-with-geometry
 * preload, and the per-interaction batch filter (which transfers the resident batch
 * rather than re-fetching the file).
 *
 * The timeout is widened because a large transcripts decode legitimately runs tens of
 * seconds in the worker; on timeout the client falls back to the main thread rather
 * than failing, so a slow worker degrades instead of breaking.
 *
 * `codec: false` because MDV starts the zarr chunk worker itself via
 * `ensureChunkWorker`, not because it wants it off.
 *
 * Returns whether the parquet worker is actually running — `false` outside a browser
 * or after a failed load — which is what UI depending on the no-fallback scan should
 * gate on. `isParquetWorkerEnabled()` from core answers the same question later.
 */
export function ensureParquetWorker(): boolean {
    // `ensureWorkers` is idempotent and attempts the parquet worker once per page, so
    // this needs no latch of its own.
    const { parquet } = ensureWorkers({
        codec: false,
        parquet: { workerUrl: parquetWorkerUrl, requestTimeoutMs: 120_000 },
    });
    return parquet;
}
