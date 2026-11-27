/**
 * PROTOTYPE (layout lab) — throwaway. The workspace prototype's own grid engine hosting real
 * charts, with the prototype's gesture handling: the card lifts and follows the pointer, a
 * placeholder shows where it will land, other cards slide out of the way, the grid grows and
 * scrolls while dragging near an edge, and charts redraw once on drop. Esc cancels. Alt+arrows
 * move the focused card and Shift+Alt+arrows resize it. Editing only at 12 columns (ADR 0026).
 */
import { useEffect, useLayoutEffect, useMemo, useRef, useState } from "react";
import ChartHost from "./ChartHost";
import {
    COLS,
    colsFor,
    compact,
    GAP_PX,
    type Item,
    type Layout12,
    layoutAt,
    ROW_PX,
    type Rect,
    toLayout,
    toPx,
} from "./engine";
import { gestureStart, gestureStop, type LabChart, type Measurement, meterStart, meterStop } from "./labShared";

type Props = {
    charts: LabChart[];
    layout: Layout12;
    onLayoutChange: (layout: Layout12) => void;
    onMeasured: (m: Measurement) => void;
};

type Box = { left: number; top: number; width: number; height: number };
type Gesture = {
    id: string;
    el: HTMLElement;
    kind: "move" | "resize";
    dir: string;
    sx: number;
    sy: number;
    px: number;
    py: number;
    scroll0: number;
    base: Item[];
    it: Item;
    box: Box;
    started: boolean;
    preview: Item[] | null;
    key: string;
    lastTick?: number;
};

const CANCEL = ".ciview-chart-icon, .ciview-menu-icon, input, select, button, textarea, a";
/** extra rows below the last card while dragging, so there is somewhere to drag to */
const DRAG_SLACK_ROWS = 6;

function place(el: HTMLElement, b: Box) {
    el.style.transform = `translate(${b.left}px, ${b.top}px)`;
    el.style.width = `${b.width}px`;
    el.style.height = `${b.height}px`;
}

export default function InHouseMode({ charts, layout, onLayoutChange, onMeasured }: Props) {
    const scrollRef = useRef<HTMLDivElement>(null);
    const gridRef = useRef<HTMLDivElement>(null);
    const phRef = useRef<HTMLDivElement>(null);
    const cards = useRef(new Map<string, HTMLDivElement>());
    const gesture = useRef<Gesture | null>(null);
    const raf = useRef(0);
    const [width, setWidth] = useState(0);

    useLayoutEffect(() => {
        const el = scrollRef.current;
        if (!el) return;
        const measure = () => setWidth(el.clientWidth);
        const ro = new ResizeObserver(measure);
        ro.observe(el);
        measure();
        return () => ro.disconnect();
    }, []);

    const cols = width ? colsFor(width) : COLS;
    const editable = cols === COLS;
    const items = useMemo(() => layoutAt(layout, cols), [layout, cols]);
    const boxOf = (r: Rect): Box => toPx(r, width, cols);
    const bottomOf = (list: Item[]) => Math.max(0, ...list.map((it) => boxOf(it).top + boxOf(it).height)) + GAP_PX;
    (window as unknown as { __labCols?: number }).__labCols = cols;

    const setGridHeight = (list: Item[], slack: number) => {
        const grid = gridRef.current;
        const scroll = scrollRef.current;
        if (!grid || !scroll) return;
        grid.style.height = `${Math.max(bottomOf(list) + slack, scroll.clientHeight)}px`;
    };

    /** Recomputes the drop target from the pointer and the scroll position, and previews it. */
    const update = () => {
        const g = gesture.current;
        const scroll = scrollRef.current;
        if (!g || !scroll) return;
        const dx = g.px - g.sx;
        const dy = g.py - g.sy + (scroll.scrollTop - g.scroll0);
        if (!g.started) {
            if (Math.hypot(dx, dy) < 4) return;
            g.started = true;
            g.el.classList.add("lift");
            gridRef.current?.classList.add("moving");
            gestureStart();
            meterStart();
            raf.current = requestAnimationFrame(autoScroll);
        }
        const colStep = (width - GAP_PX * (cols + 1)) / cols + GAP_PX;
        const rowStep = ROW_PX + GAP_PX;
        const t: Item = { ...g.it };
        if (g.kind === "move") {
            place(g.el, { ...g.box, left: g.box.left + dx, top: g.box.top + dy });
            t.x = Math.max(0, Math.min(cols - t.w, Math.round((g.box.left + dx - GAP_PX) / colStep)));
            t.y = Math.max(0, Math.round((g.box.top + dy - GAP_PX) / rowStep));
        } else {
            const west = g.dir.includes("w");
            const dW = west ? -dx : g.dir.includes("e") ? dx : 0;
            const dH = g.dir.includes("s") ? dy : 0;
            const wPx = Math.max(80, g.box.width + dW);
            const hPx = Math.max(80, g.box.height + dH);
            place(g.el, { left: west ? g.box.left + g.box.width - wPx : g.box.left, top: g.box.top, width: wPx, height: hPx });
            const w = Math.max(2, Math.round((wPx + GAP_PX) / colStep));
            t.h = Math.max(3, Math.round((hPx + GAP_PX) / rowStep));
            if (west) {
                const right = g.it.x + g.it.w;
                t.w = Math.min(w, right);
                t.x = right - t.w;
            } else {
                t.w = Math.min(w, cols - g.it.x);
            }
        }
        const key = `${t.x},${t.y},${t.w},${t.h}`;
        if (key === g.key) return;
        g.key = key;
        const pushed = compact(
            g.base.map((o) => (o.i === g.id ? t : o)),
            cols,
            new Set([g.id]),
        );
        const settled = compact(pushed, cols);
        g.preview = settled;
        for (const o of settled) {
            if (o.i === g.id) {
                if (phRef.current) place(phRef.current, boxOf(o));
                continue;
            }
            const el = cards.current.get(o.i);
            if (el) place(el, boxOf(o));
        }
        setGridHeight(settled, DRAG_SLACK_ROWS * rowStep);
    };

    /** Scrolls while the pointer is held near the top or bottom edge, faster the closer it is. */
    const autoScroll = () => {
        const g = gesture.current;
        const scroll = scrollRef.current;
        if (!g || !scroll) return;
        const r = scroll.getBoundingClientRect();
        const zone = 56;
        // up to ~600 px a second, whatever the frame rate
        const now = performance.now();
        const dt = Math.min(50, now - (g.lastTick ?? now));
        g.lastTick = now;
        const speed = 0.6 * dt;
        let v = 0;
        if (g.py > r.bottom - zone) v = Math.min(1, (g.py - (r.bottom - zone)) / zone) * speed;
        else if (g.py < r.top + zone) v = -Math.min(1, (r.top + zone - g.py) / zone) * speed;
        if (v) {
            const before = scroll.scrollTop;
            scroll.scrollTop += v;
            if (scroll.scrollTop !== before) update();
        }
        raf.current = requestAnimationFrame(autoScroll);
    };

    const onMove = (e: PointerEvent) => {
        const g = gesture.current;
        if (!g) return;
        g.px = e.clientX;
        g.py = e.clientY;
        update();
    };

    const end = (commit: boolean) => {
        window.removeEventListener("pointermove", onMove);
        window.removeEventListener("pointerup", onUp);
        window.removeEventListener("pointercancel", onCancel);
        window.removeEventListener("keydown", onKey);
        cancelAnimationFrame(raf.current);
        const g = gesture.current;
        gesture.current = null;
        if (!g) return;
        g.el.classList.remove("lift");
        gridRef.current?.classList.remove("moving");
        if (!g.started) return;
        const final = commit && g.preview ? g.preview : g.base;
        // React will not touch a card whose committed place equals what it last rendered, so put
        // every card where it belongs now, rather than leaving the dragged one where it was let go
        for (const o of final) {
            const el = cards.current.get(o.i);
            if (el) place(el, boxOf(o));
        }
        setGridHeight(final, 0);
        if (commit && g.preview) onLayoutChange(toLayout(g.preview));
        gestureStop();
        void meterStop(`in-house ${g.kind}${commit ? "" : " (cancelled)"} @ ${cols} cols`).then(onMeasured);
    };
    const onUp = () => end(true);
    const onCancel = () => end(false);
    const onKey = (e: KeyboardEvent) => {
        if (e.key === "Escape") end(false);
    };

    const begin = (e: React.PointerEvent, id: string, kind: "move" | "resize", dir = "") => {
        if (!editable || e.button !== 0) return;
        const it = items.find((o) => o.i === id);
        const el = cards.current.get(id);
        const scroll = scrollRef.current;
        if (!it || !el || !scroll) return;
        e.preventDefault();
        gesture.current = {
            id,
            el,
            kind,
            dir,
            sx: e.clientX,
            sy: e.clientY,
            px: e.clientX,
            py: e.clientY,
            scroll0: scroll.scrollTop,
            base: items,
            it,
            box: boxOf(it),
            started: false,
            preview: null,
            key: "",
        };
        window.addEventListener("pointermove", onMove);
        window.addEventListener("pointerup", onUp);
        window.addEventListener("pointercancel", onCancel);
        window.addEventListener("keydown", onKey);
    };

    // a gesture still running when the engine unmounts must not leave redraws held back
    useEffect(() => () => end(false), []);

    const nudge = (id: string, dx: number, dy: number, dw: number, dh: number) => {
        const it = items.find((o) => o.i === id);
        if (!it || !editable) return;
        const t = {
            ...it,
            x: Math.max(0, Math.min(cols - it.w, it.x + dx)),
            y: Math.max(0, it.y + dy),
            w: Math.max(2, Math.min(cols - it.x, it.w + dw)),
            h: Math.max(3, it.h + dh),
        };
        const next = compact(
            compact(
                items.map((o) => (o.i === id ? t : o)),
                cols,
                new Set([id]),
            ),
            cols,
        );
        onLayoutChange(toLayout(next));
    };

    const onCardKey = (e: React.KeyboardEvent, id: string) => {
        if (!e.altKey || e.target !== e.currentTarget) return;
        const d = { ArrowLeft: [-1, 0], ArrowRight: [1, 0], ArrowUp: [0, -1], ArrowDown: [0, 1] }[e.key];
        if (!d) return;
        e.preventDefault();
        if (e.shiftKey) nudge(id, 0, 0, d[0], d[1]);
        else nudge(id, d[0], d[1], 0, 0);
    };

    const height = width ? bottomOf(items) : 0;
    return (
        <div ref={scrollRef} className="lab-scroll lab-canvas">
            <div ref={gridRef} className="lab-grid" style={{ height, minHeight: "100%" }} data-cols={cols}>
                <div ref={phRef} className="lab-ph" />
                {width > 0 &&
                    charts.map((c) => {
                        const it = items.find((o) => o.i === c.id);
                        if (!it) return null;
                        const b = boxOf(it);
                        return (
                            <div
                                key={c.id}
                                ref={(el) => {
                                    if (el) cards.current.set(c.id, el);
                                    else cards.current.delete(c.id);
                                }}
                                className="lab-card"
                                style={{
                                    transform: `translate(${b.left}px, ${b.top}px)`,
                                    width: b.width,
                                    height: b.height,
                                }}
                                data-id={c.id}
                                tabIndex={0}
                                aria-label={`${c.title}: column ${it.x + 1}, row ${it.y + 1}, ${it.w} wide, ${it.h} tall`}
                                onKeyDown={(e) => onCardKey(e, c.id)}
                                onPointerDown={(e) => {
                                    const t = e.target as HTMLElement;
                                    if (t.closest(".ciview-chart-titlebar") && !t.closest(CANCEL)) begin(e, c.id, "move");
                                }}
                            >
                                <ChartHost chart={c.chart} />
                                {editable &&
                                    ["e", "s", "se", "w", "sw"].map((d) => (
                                        <div
                                            key={d}
                                            className={`lab-rz lab-rz-${d}`}
                                            onPointerDown={(e) => {
                                                e.stopPropagation();
                                                begin(e, c.id, "resize", d);
                                            }}
                                        />
                                    ))}
                            </div>
                        );
                    })}
            </div>
        </div>
    );
}
