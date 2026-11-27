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
 * | `association.rows` | **absent** during a cosmetic edit or filter; once per geometry or key-column load | Feature ids are being re-resolved to DataStore rows (`table_rows.ts`). It depends only on geometry, table association and key-column loads — never on layer props. |
 * | `adapter.*` | one per viewer render | These are cheap; the `count` is the signal, and it tracks how often the adapter re-ran. |
 * | `render:spatial.viewer` / `.canvas` | one per interaction step, under ~16ms | More than one per step means something downstream is setting state during the render that follows an edit. |
 * | `loaded:first` | one per chart opened | Time from opening the chart to the first frame with every layer fully loaded. Long means slow metadata, chunk fetches or decoding; check `raster_health.ts` issues for the image. |
 * | `loaded:settle` | one per pan or zoom that needed new tiles | Time from the last view change to every visible tile loaded. Long with a high count usually means oversized chunks. |
 */
import { recordBrowserTiming } from "@/utils/browserTelemetry";

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
        const table: Record<string, { count: number; avgMs: number; maxMs: number; totalMs: number }> = {};
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

/**
 * How long the canvas took to draw every layer fully loaded, from when the chart opened
 * (`first`) or from the last pan or zoom that needed new tiles (`settle`). See
 * `load_timing.ts` for how the end is found.
 *
 * Unlike the labels above, these also go to browser telemetry when the deployment has it
 * on (ADR 0025), because they are what people wait for, rather than what the code does.
 */
export function recordSpatialLoaded(kind: "first" | "settle", ms: number, context: Record<string, string>) {
    if (ENABLED) bump(`loaded:${kind}`, ms);
    recordBrowserTiming(`mdv_spatial_${kind}_load`, { duration_ms: ms }, context);
}

/** React `<Profiler onRender>` adapter — records committed render duration per id. */
export function onSpatialProfilerRender(id: string, _phase: unknown, actualDuration: number) {
    if (!ENABLED) return;
    bump(`render:${id}`, actualDuration);
}

/**
 * A deck to time whole frames of, as `__benchSpatialDeck(label)`. Also records each frame's
 * CPU time under `<label>.draw`. `pan(step)` moves the view `step` device pixels from where
 * it is (0 puts it back), so a frame can be forced to draw from scratch. No-op when disabled.
 *
 * ## Frame benchmarks
 *
 * `await __benchSpatialDeck("densityGrid", { frames: 60, pan: true })` draws `frames` frames
 * synchronously, each waiting for the GPU before the next, and returns ms per frame. With
 * `pan`, every frame moves the view a pixel, so cached drawing (an isolated group's target)
 * is redrawn once per frame; without, nothing changed, so only what cannot be cached draws.
 * Labels recorded meanwhile, such as `densityGrid.targetRender` (one per isolated-group
 * target drawn), say how much was drawn.
 */
type PerfDeck = {
    redraw: (reason?: string) => void;
    setProps: (props: { onBeforeRender?: () => void; onAfterRender?: () => void }) => void;
    getCanvasContext: () => { device: { handle: unknown; submit: () => void } } | null;
};

type PerfDeckEntry = { deck: PerfDeck; pan: (step: number) => void };
const perfDecks = new Map<string, PerfDeckEntry>();

function hasGPUQueue(handle: unknown): handle is { queue: { onSubmittedWorkDone: () => Promise<void> } } {
    return (
        typeof handle === "object" &&
        handle !== null &&
        "queue" in handle &&
        typeof handle.queue === "object" &&
        handle.queue !== null &&
        "onSubmittedWorkDone" in handle.queue
    );
}

/** Waits for the GPU to finish the work submitted so far; says how, as a check on the timings. */
async function waitForGPU(deck: PerfDeck) {
    const device = deck.getCanvasContext()?.device;
    if (!device) return "no device";
    device.submit();
    const { handle } = device;
    if (hasGPUQueue(handle)) {
        await handle.queue.onSubmittedWorkDone();
        return "webgpu";
    }
    if (typeof WebGL2RenderingContext !== "undefined" && handle instanceof WebGL2RenderingContext) {
        // finish() returns early in some browsers; reading a pixel cannot
        handle.readPixels(0, 0, 1, 1, handle.RGBA, handle.UNSIGNED_BYTE, new Uint8Array(4));
        return "webgl";
    }
    return "none";
}

async function benchSpatialDeck(label: string, { frames = 60, pan = true } = {}) {
    const entry = perfDecks.get(label);
    if (!entry) throw new Error(`no deck registered as "${label}"; have ${[...perfDecks.keys()].join(", ")}`);
    const { deck } = entry;
    const sync = await waitForGPU(deck);
    const t0 = performance.now();
    for (let i = 0; i < frames; i++) {
        if (pan) entry.pan(i % 2 === 0 ? 1 : 0);
        deck.redraw("perf bench");
        await waitForGPU(deck);
    }
    const msPerFrame = (performance.now() - t0) / frames;
    entry.pan(0);
    deck.redraw("perf bench");
    return { msPerFrame: +msPerFrame.toFixed(2), sync };
}

export function registerSpatialPerfDeck(label: string, deck: PerfDeck, pan: (step: number) => void) {
    if (!ENABLED) return () => {};
    let frameStart = 0;
    deck.setProps({
        onBeforeRender: () => {
            frameStart = performance.now();
        },
        onAfterRender: () => bump(`${label}.draw`, performance.now() - frameStart),
    });
    perfDecks.set(label, { deck, pan });
    (window as unknown as { __benchSpatialDeck?: typeof benchSpatialDeck }).__benchSpatialDeck = benchSpatialDeck;
    return () => {
        if (perfDecks.get(label)?.deck === deck) perfDecks.delete(label);
    };
}
