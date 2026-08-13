import { enablePointsWorker, setPointsWorkerRequestTimeout } from "@spatialdata/core";

// Vite bundles the core points worker and hands back a runtime URL. The explicit URL
// is required: `enablePointsWorker()` with no options builds
// `new URL("./points-worker.js", import.meta.url)` behind a `@vite-ignore`, which
// resolves relative to wherever the dep was served from (`node_modules/.vite/deps` in
// dev) and 404s.
import pointsWorkerUrl from "@spatialdata/core/points-worker?worker&url";

let pointsWorkerEnabled = false;

/**
 * Enable the core points worker before any spatial points work.
 *
 * Not optional for anything past the first capped read. Core's `defaultEnabled` is
 * `false`, and `loadPointsMatchingFeatureCodes` — the feature-index scan that fetches a
 * selected feature's points from beyond the resident window — throws outright without
 * it. Until @spatialdata/core 0.8.0 this could not be done at all: the published
 * worker entry was a CommonJS file in an ESM package, so the Worker died on
 * `require is not defined` (SpatialData.js#148). That is why the pin floor is 0.8.0.
 *
 * It also moves the heavy decodes off the main thread: the codes-with-geometry
 * preload, and the per-interaction batch filter (which transfers the resident batch
 * rather than re-fetching the file).
 *
 * The timeout is widened because a large transcripts decode legitimately runs tens of
 * seconds in the worker; on timeout the client falls back to the main thread rather
 * than failing, so a slow worker degrades instead of breaking.
 */
export function ensurePointsWorker() {
    if (pointsWorkerEnabled) return;
    enablePointsWorker({ workerUrl: pointsWorkerUrl });
    setPointsWorkerRequestTimeout(120_000);
    pointsWorkerEnabled = true;
}
