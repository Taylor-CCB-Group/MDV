/**
 * Opt-in perf instrumentation for the spatial image render path.
 *
 * Off by default and zero-cost when off. Enable from the browser console:
 *
 *   localStorage.MDV_SPATIAL_PERF = "1"   // then reload
 *
 * Then interact (e.g. drag a contrast slider for ~2s) and read the breakdown on
 * `window.__spatialPerf`, or watch the throttled `console.table`. `count` is how
 * many times each measured step ran during the capture — for adapter steps that
 * equals the number of `SpatialDataViewer` re-renders, which is the key signal:
 * a cosmetic image edit should not re-render the whole viewer on every frame.
 *
 * Call `window.__resetSpatialPerf()` between captures.
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
