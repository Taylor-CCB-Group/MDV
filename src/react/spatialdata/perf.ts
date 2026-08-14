/**
 * Opt-in perf instrumentation for the spatial render path — the single place this
 * exists. If you are adding a timer somewhere in `src/react/spatialdata`, add a
 * label here rather than a second mechanism.
 *
 * Off by default and zero-cost when off. Enable from the browser console:
 *
 *   localStorage.MDV_SPATIAL_PERF = "1"   // then reload
 *   delete localStorage.MDV_SPATIAL_PERF  // then reload
 *
 * **The flag is sticky.** It lives in localStorage, so it survives reloads, branch
 * switches and days off, and `ENABLED` is read once at module load — toggling it
 * needs a reload either way. A capture left on from last week is the classic way to
 * measure the wrong thing, or to chase a "slowdown" that is the instrumentation. So
 * when it is on it says so, loudly, on every load. Do not silence that banner.
 *
 * Interact (drag a slider for ~2s), then read `window.__spatialPerf` or the
 * throttled `console.table`. `window.__resetSpatialPerf()` clears between captures.
 *
 * ## Reading a capture
 *
 * `count` is usually the more informative number: it says how many times a step ran
 * during the interaction, and most of the regressions this has caught were "ran at
 * all" rather than "ran slowly".
 *
 * | Label | Healthy | What it means when it is not |
 * |---|---|---|
 * | `shapes.loadRenderData` | **absent**, or one per element per session | Geometry is being re-decoded. Parquet read + WKB decode of every polygon, ~0.5s each. Nothing cosmetic should ever produce one — see `useElementKeys` in `table_association.ts`. |
 * | `association.project` | high `count`, `avgMs` ≈ 0 | High `avgMs` means the per-feature pass is rebuilding rather than reusing its cache. Expected on a colour or filter change; on an opacity drag it is a bug. |
 * | `adapter.*` | one per viewer render | These are cheap; the `count` is the signal, and it tracks how often the adapter re-ran. |
 * | `render:spatial.viewer` / `.canvas` | one per interaction step, under ~16ms | More than one per step means something downstream is setting state during the render that follows an edit. |
 */

type Stat = { count: number; total: number; max: number };

const ENABLED =
    typeof window !== "undefined" &&
    (() => {
        try {
            return window.localStorage?.getItem("MDV_SPATIAL_PERF") === "1";
        } catch {
            return false;
        }
    })();

const stats = new Map<string, Stat>();
let flushTimer: ReturnType<typeof setTimeout> | null = null;

function bump(label: string, ms: number) {
    let s = stats.get(label);
    if (!s) {
        s = { count: 0, total: 0, max: 0 };
        stats.set(label, s);
    }
    s.count++;
    s.total += ms;
    if (ms > s.max) s.max = ms;
    scheduleFlush();
}

function scheduleFlush() {
    if (flushTimer) return;
    flushTimer = setTimeout(() => {
        flushTimer = null;
        const table: Record<
            string,
            { count: number; avgMs: number; maxMs: number; totalMs: number }
        > = {};
        for (const [label, s] of stats) {
            table[label] = {
                count: s.count,
                avgMs: +(s.total / Math.max(1, s.count)).toFixed(2),
                maxMs: +s.max.toFixed(2),
                totalMs: +s.total.toFixed(2),
            };
        }
        (window as unknown as { __spatialPerf?: unknown }).__spatialPerf = table;
        console.table(table);
    }, 500);
}

if (ENABLED && typeof window !== "undefined") {
    (window as unknown as { __resetSpatialPerf?: () => void }).__resetSpatialPerf = () => {
        stats.clear();
        (window as unknown as { __spatialPerf?: unknown }).__spatialPerf = {};
    };
    (window as unknown as { __disableSpatialPerf?: () => void }).__disableSpatialPerf = () => {
        window.localStorage?.removeItem("MDV_SPATIAL_PERF");
        console.info("[MDV spatial perf] disabled — reload to take effect.");
    };
    // The flag is sticky and there is no UI for it, so the only thing standing
    // between a stale capture and someone believing it is this banner. It is
    // deliberately noisy: whoever set this may have been a different person, on a
    // different branch, weeks ago.
    console.warn(
        "%c[MDV spatial perf] INSTRUMENTATION IS ON%c\n" +
            "localStorage.MDV_SPATIAL_PERF is set, so the spatial render path is being timed.\n" +
            "Timings include this overhead, and any perf comparison against a machine without\n" +
            "the flag is not like-for-like. To turn it off:  __disableSpatialPerf()  then reload.\n" +
            "Labels and what a healthy capture looks like: src/react/spatialdata/perf.ts",
        "background:#b45309;color:#fff;font-weight:bold;padding:2px 6px;border-radius:3px",
        "",
    );
}

export const spatialPerfEnabled = ENABLED;

/** Record a pre-measured duration (ms) under `label`. No-op when disabled. */
export function recordSpatialPerf(label: string, ms: number) {
    if (!ENABLED) return;
    bump(label, ms);
}

/** Time `fn` under `label`. Returns `fn`'s result. No-op wrapper when disabled. */
export function measureSpatial<T>(label: string, fn: () => T): T {
    if (!ENABLED) return fn();
    const t0 = performance.now();
    try {
        return fn();
    } finally {
        bump(label, performance.now() - t0);
    }
}

/** React `<Profiler onRender>` adapter — records committed render duration per id. */
export function onSpatialProfilerRender(
    id: string,
    _phase: unknown,
    actualDuration: number,
) {
    if (!ENABLED) return;
    bump(`render:${id}`, actualDuration);
}
