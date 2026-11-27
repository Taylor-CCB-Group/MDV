/**
 * PROTOTYPE (layout lab) — throwaway. The grid model both engines are judged against
 * (ADR 0026): one saved 12-column layout in grid units; 6- and 1-column layouts derived
 * from it in reading order and never stored. The in-house engine also uses `compact`.
 */

export type Rect = { x: number; y: number; w: number; h: number };
export type Item = Rect & { i: string };
/** The saved shape: chart id to its place on the 12-column grid. */
export type Layout12 = Record<string, Rect>;

export const COLS = 12;
export const ROW_PX = 50;
export const GAP_PX = 8;

/** 12 columns from 960 px, 6 from 600 px, 1 below (ADR 0026). */
export const colsFor = (width: number) => (width >= 960 ? 12 : width >= 600 ? 6 : 1);

const clamp = (v: number, lo: number, hi: number) => Math.max(lo, Math.min(hi, v));

export const toItems = (layout: Layout12): Item[] => Object.entries(layout).map(([i, r]) => ({ i, ...r }));
export const toLayout = (items: Item[]): Layout12 =>
    Object.fromEntries(items.map(({ i, x, y, w, h }) => [i, { x, y, w, h }]));

const byReadingOrder = (a: Item, b: Item) => a.y - b.y || a.x - b.x || (a.i < b.i ? -1 : 1);
export const hit = (a: Item, b: Item) =>
    a.i !== b.i && a.x < b.x + b.w && b.x < a.x + a.w && a.y < b.y + b.h && b.y < a.y + a.h;
const hitAny = (a: Item, list: Item[]) => list.some((b) => hit(a, b));

/** Pulls items up into free space. Pinned items (the one being dragged) stay put; others flow around them. */
export function compact(items: Item[], cols: number, pinned: Set<string> = new Set()): Item[] {
    const out = items.map((o) => ({ ...o }));
    const placed = out.filter((o) => pinned.has(o.i));
    for (const it of out.filter((o) => !pinned.has(o.i)).sort(byReadingOrder)) {
        it.w = Math.min(it.w, cols);
        it.x = clamp(it.x, 0, cols - it.w);
        let y = Math.max(0, it.y);
        while (y > 0 && !hitAny({ ...it, y: y - 1 }, placed)) y--;
        while (hitAny({ ...it, y }, placed)) y++;
        it.y = y;
        placed.push(it);
    }
    return out;
}

/** A narrower layout: the 12-column one in reading order, widths scaled, rows wrapped. */
export function flow(layout: Layout12, cols: number): Item[] {
    const out: Item[] = [];
    let cx = 0;
    let cy = 0;
    let rowH = 0;
    for (const it of toItems(layout).sort(byReadingOrder)) {
        const w = cols === 1 ? 1 : clamp(Math.round((it.w * cols) / COLS), Math.min(3, cols), cols);
        if (cx + w > cols) {
            cx = 0;
            cy += rowH;
            rowH = 0;
        }
        out.push({ i: it.i, x: cx, y: cy, w, h: it.h });
        cx += w;
        rowH = Math.max(rowH, it.h);
    }
    return out;
}

/** The layout shown at `cols` columns. */
export const layoutAt = (layout: Layout12, cols: number): Item[] =>
    cols === COLS ? compact(toItems(layout), COLS) : flow(layout, cols);

/** Pixel box of an item, for a grid `width` px wide. */
export function toPx(item: Rect, width: number, cols: number) {
    const colW = (width - GAP_PX * (cols + 1)) / cols;
    return {
        left: GAP_PX + item.x * (colW + GAP_PX),
        top: GAP_PX + item.y * (ROW_PX + GAP_PX),
        width: item.w * colW + (item.w - 1) * GAP_PX,
        height: item.h * ROW_PX + (item.h - 1) * GAP_PX,
        colW,
    };
}

/** Same layout, ignoring key order, so a reload can be checked for drift. */
export function sameLayout(a: Layout12, b: Layout12) {
    const ka = Object.keys(a).sort();
    const kb = Object.keys(b).sort();
    if (ka.join() !== kb.join()) return false;
    return ka.every((k) => ["x", "y", "w", "h"].every((p) => a[k][p as keyof Rect] === b[k][p as keyof Rect]));
}
