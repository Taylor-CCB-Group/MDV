/**
 * PROTOTYPE (layout lab) — throwaway. Shared pieces for the layout engine bake-off
 * (gridstack 14 vs the in-house engine) on a project's real charts.
 */
import type BaseChart from "@/charts/BaseChart";
import type ChartManager from "@/charts/ChartManager";
import { COLS, compact, GAP_PX, type Item, type Layout12, ROW_PX, toLayout } from "./engine";

export type LabChart = {
    id: string;
    title: string;
    dsName: string;
    chart: BaseChart<any>;
};

export function collectCharts(cm: ChartManager): LabChart[] {
    return Object.entries(cm.charts).map(([id, state]) => ({
        id,
        title: state.chart.config.title || state.chart.config.type,
        dsName: state.dataSource.name,
        chart: state.chart,
    }));
}

/** gridstack rows in MDV are ~150 px (`GridstackManager.getGrid`), about three lab rows. */
const GS_ROW_PX = 150;

/**
 * Converts a pane view as it opens (ADR 0026): pane by pane in the old left-to-right order,
 * within a pane top to bottom then left to right; each chart keeps its width relative to its
 * pane and is packed row by row into 12 columns, then pulled up so no gaps are left.
 */
export function convertPaneView(cm: ChartManager, charts: LabChart[]): Layout12 {
    const panes = (cm.viewData?.dataSources ?? {}) as Record<string, { layout?: string } | undefined>;
    const order = [...Object.keys(panes)];
    for (const c of charts) if (!order.includes(c.dsName)) order.push(c.dsName);

    type Rel = { id: string; top: number; left: number; w: number; hPx: number };
    const sequence: Rel[] = [];
    for (const ds of order) {
        const inPane = charts.filter((c) => c.dsName === ds);
        const grid = panes[ds]?.layout === "gridstack";
        // pixel layouts were saved at whatever width the author's screen had: measure from that
        const paneW = Math.max(
            600,
            ...inPane.map(({ chart }) => (chart.config.position?.[0] ?? 0) + (chart.config.size?.[0] ?? 0)),
        );
        const rels = inPane.map(({ id, chart }): Rel => {
            const cfg = chart.config;
            if (grid && cfg.gsposition && cfg.gssize) {
                const [gx, gy] = cfg.gsposition;
                const [gw, gh] = cfg.gssize;
                return { id, left: gx / COLS, top: gy * GS_ROW_PX, w: gw / COLS, hPx: gh * GS_ROW_PX };
            }
            const [left, top] = cfg.position ?? [0, 0];
            const [width, height] = cfg.size ?? [400, 300];
            return { id, left: left / paneW, top, w: width / paneW, hPx: height };
        });
        rels.sort((a, b) => a.top - b.top || a.left - b.left || (a.id < b.id ? -1 : 1));
        sequence.push(...rels);
    }

    const items: Item[] = [];
    let cx = 0;
    let cy = 0;
    let rowH = 0;
    for (const r of sequence) {
        const w = Math.max(2, Math.min(COLS, Math.round(r.w * COLS)));
        const h = Math.max(3, Math.round((r.hPx + GAP_PX) / (ROW_PX + GAP_PX)));
        if (cx + w > COLS) {
            cx = 0;
            cy += rowH;
            rowH = 0;
        }
        items.push({ i: r.id, x: cx, y: cy, w, h });
        cx += w;
        rowH = Math.max(rowH, h);
    }
    return toLayout(compact(items, COLS));
}

/* ---------------- saving (stands in for views_v2.json) ---------------- */

const storageKey = (cm: ChartManager) =>
    `mdv.layoutLab.v2:${window.location.pathname}:${cm.viewManager?.current_view ?? ""}`;

export function loadSaved(cm: ChartManager): Layout12 | null {
    try {
        const raw = window.localStorage.getItem(storageKey(cm));
        const parsed = raw ? JSON.parse(raw) : null;
        return parsed?.version === 2 ? parsed.layout : null;
    } catch {
        return null;
    }
}

export function save(cm: ChartManager, layout: Layout12) {
    try {
        window.localStorage.setItem(storageKey(cm), JSON.stringify({ version: 2, layout }));
    } catch (error) {
        console.warn("layout lab: could not save", error);
    }
}

export function clearSaved(cm: ChartManager) {
    try {
        window.localStorage.removeItem(storageKey(cm));
    } catch {
        /* nothing to clear */
    }
}

/* ---------------- gestures ---------------- */

/**
 * While a chart is being dragged or resized, hosts hold back chart redraws and do them once
 * when the gesture ends, as the prototype does. Every engine reports its gestures here.
 */
let gestureActive = false;
const gestureEnd = new Set<() => void>();
export const isGestureActive = () => gestureActive;
export function gestureStart() {
    gestureActive = true;
}
export function gestureStop() {
    gestureActive = false;
    for (const fn of [...gestureEnd]) fn();
}
export function onGestureEnd(fn: () => void) {
    gestureEnd.add(fn);
    return () => gestureEnd.delete(fn);
}

/* ---------------- interaction meter ---------------- */

export type Measurement = {
    label: string;
    frames: number;
    avgFps: number;
    worstFrameMs: number;
    longFrames: number;
    setSizeCalls: number;
};

let setSizeCalls = 0;
export function countSetSize() {
    setSizeCalls++;
}

let raf = 0;
let last = 0;
let deltas: number[] = [];
let startCalls = 0;

export function meterStart() {
    cancelAnimationFrame(raf);
    deltas = [];
    startCalls = setSizeCalls;
    last = performance.now();
    const tick = (t: number) => {
        deltas.push(t - last);
        last = t;
        raf = requestAnimationFrame(tick);
    };
    raf = requestAnimationFrame(tick);
}

/** Stops the meter; `settleMs` keeps counting resizes the debounced hosts fire after release. */
export function meterStop(label: string, settleMs = 250): Promise<Measurement> {
    cancelAnimationFrame(raf);
    const frameDeltas = deltas.slice(1);
    return new Promise((resolve) => {
        setTimeout(() => {
            const total = frameDeltas.reduce((a, b) => a + b, 0);
            const m: Measurement = {
                label,
                frames: frameDeltas.length,
                avgFps: frameDeltas.length ? Math.round((1000 * frameDeltas.length) / total) : 0,
                worstFrameMs: Math.round(Math.max(0, ...frameDeltas)),
                longFrames: frameDeltas.filter((d) => d > 50).length,
                setSizeCalls: setSizeCalls - startCalls,
            };
            const w = window as unknown as { __layoutLab?: Measurement[] };
            w.__layoutLab = [...(w.__layoutLab ?? []), m];
            console.info("[layout lab]", m);
            resolve(m);
        }, settleMs);
    });
}
