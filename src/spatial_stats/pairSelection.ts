/** A chosen pair as stored in config: category labels, anchor A → neighbour B (g is asymmetric). */
export type LabelPair = { a: string; b: string };

const same = (p: LabelPair, q: LabelPair) => p.a === q.a && p.b === q.b;
export const hasPair = (pairs: readonly LabelPair[], p: LabelPair) => pairs.some((q) => same(p, q));

/** Configured pairs as type ids of this result; pairs naming a category it doesn't have are dropped,
 *  and the rest keep their position in config (which their colour follows). */
export function resolvePairs(labels: readonly string[], pairs: readonly LabelPair[]) {
    const out: { a: number; b: number; index: number }[] = [];
    pairs.forEach((p, index) => {
        const a = labels.indexOf(p.a);
        const b = labels.indexOf(p.b);
        if (a >= 0 && b >= 0) out.push({ a, b, index });
    });
    return out;
}

export function addPair(pairs: readonly LabelPair[], p: LabelPair): LabelPair[] {
    return hasPair(pairs, p) ? [...pairs] : [...pairs, p];
}

export function removePair(pairs: readonly LabelPair[], p: LabelPair): LabelPair[] {
    return pairs.filter((q) => !same(p, q));
}

/** `pairs` with A → B added, or removed if already there. */
export function togglePair(pairs: readonly LabelPair[], a: string, b: string): LabelPair[] {
    const p = { a, b };
    return hasPair(pairs, p) ? removePair(pairs, p) : addPair(pairs, p);
}

/**
 * Brushing pairs: a pointer-down starts a stroke, every further cell the pointer passes over continues
 * it. A plain press starts afresh with just that pair and the stroke adds. With a modifier the stroke
 * edits the current pairs: it removes if the first cell was already chosen, otherwise adds — so one
 * stroke never flickers a pair on and off. A click is a stroke of one cell.
 */
export type PaintMode = "add" | "remove";

export function startPaint(
    pairs: readonly LabelPair[],
    p: LabelPair,
    modifier: boolean,
): { pairs: LabelPair[]; mode: PaintMode } {
    if (!modifier) return { pairs: [p], mode: "add" };
    return hasPair(pairs, p)
        ? { pairs: removePair(pairs, p), mode: "remove" }
        : { pairs: addPair(pairs, p), mode: "add" };
}

export function continuePaint(pairs: readonly LabelPair[], p: LabelPair, mode: PaintMode): LabelPair[] {
    return mode === "add" ? addPair(pairs, p) : removePair(pairs, p);
}

export type Cell = { i: number; j: number };

/**
 * The grid cells on the way from `from` to `to`, excluding `from`, in order — so a fast drag whose
 * pointer events are several cells apart still brushes everything it crossed.
 */
export function cellsAlong(from: Cell, to: Cell): Cell[] {
    const di = to.i - from.i;
    const dj = to.j - from.j;
    const steps = Math.max(Math.abs(di), Math.abs(dj)) * 2;
    const out: Cell[] = [];
    let last = from;
    for (let s = 1; s <= steps; s++) {
        const c = { i: Math.round(from.i + (di * s) / steps), j: Math.round(from.j + (dj * s) / steps) };
        if (c.i !== last.i || c.j !== last.j) {
            out.push(c);
            last = c;
        }
    }
    return out;
}
