import type { PcfMatrixBinnedResult } from "intraspatial/spatial/pcf";
import { type PointerEvent, useMemo, useRef } from "react";
import { type Cell, cellsAlong } from "./pairSelection";

/** `g` at one bin, or averaged over all bins. */
export type RSelection = { kind: "mean" } | { kind: "bin"; k: number };

export type MatrixView = {
    /** Type ids with cells, in matrix order. */
    types: number[];
    /** `value[i * n + j]` = g(types[i] → types[j]) at the chosen r. */
    value: Float64Array;
    pairs: Float64Array;
};

/** Reduce the K×K×B result to the present types at one r (or the r-mean). Pairs are summed over the same bins. */
export function matrixAt(res: PcfMatrixBinnedResult, r: RSelection, minCount = 1): MatrixView {
    const { nTypes: K, nBins: B, g, pairs: P, counts } = res;
    const types = [...Array(K).keys()].filter((t) => counts[t] >= minCount);
    const n = types.length;
    const value = new Float64Array(n * n);
    const pairs = new Float64Array(n * n);
    const [k0, k1] = r.kind === "mean" ? [0, B] : [r.k, r.k + 1];
    for (let i = 0; i < n; i++)
        for (let j = 0; j < n; j++) {
            const base = types[i] * K * B + types[j] * B;
            let s = 0;
            let p = 0;
            for (let k = k0; k < k1; k++) {
                s += g[base + k];
                p += P[base + k];
            }
            value[i * n + j] = s / (k1 - k0);
            pairs[i * n + j] = p;
        }
    return { types, value, pairs };
}

/** g_{a→b}(r) over all bins. */
export function pairCurve(res: PcfMatrixBinnedResult, a: number, b: number): number[] {
    const { nTypes: K, nBins: B, g } = res;
    const base = a * K * B + b * B;
    return Array.from(g.subarray(base, base + B));
}

/** Diverging colour on log2 g: blue = avoidance (g < 1), white = independence, red = attraction. */
export function gColour(g: number, maxLog2 = 2): string {
    if (!Number.isFinite(g) || g <= 0) return "rgb(40,40,40)";
    const t = Math.max(-1, Math.min(1, Math.log2(g) / maxLog2));
    const a = Math.abs(t);
    const mix = (c: number) => Math.round(255 + (c - 255) * a);
    return t >= 0 ? `rgb(${mix(200)},${mix(30)},${mix(40)})` : `rgb(${mix(30)},${mix(80)},${mix(200)})`;
}

type Pair = { a: number; b: number };

/** A brush stroke over the matrix, in type ids: `start` on pointer-down (with whether a modifier was
 *  held), then `cell` for every further cell the pointer passes over. */
export type PaintEvent = { kind: "start"; pair: Pair; modifier: boolean } | { kind: "cell"; pair: Pair };

type Props = {
    result: PcfMatrixBinnedResult;
    labels: string[];
    r: RSelection;
    /** Side of one matrix cell, px — the same for every matrix in the grid. */
    cell: number;
    /** Gutter for the type labels, and whether this matrix shows them: in a grid only the first column
     *  carries row labels and only the first row column labels, all aligned to the same cells. */
    labelSize: number;
    rowLabels: boolean;
    colLabels: boolean;
    /** Chosen pairs as type ids, with their curve colours — owned by the caller (chart config). */
    selected: readonly (Pair & { colour: string })[];
    /** The pair highlighted across the chart (hovered here, or on a curve or legend item). */
    focus: Pair | null;
    onFocus?: (pair: Pair | null) => void;
    onPaint?: (e: PaintEvent) => void;
};

/**
 * K×K heatmap of cross-PCF at one r (row = anchor type A, column = neighbour type B).
 *
 * Dragging paints pairs: a plain press starts afresh with the pressed pair and adds every cell the
 * pointer then passes over; Shift / Ctrl / ⌘ edit the current pairs instead (see `startPaint`). Chosen
 * pairs carry a thin ring in their curve colour, inside the cell so neighbours never overlap; the
 * focused pair gets a stronger one, and its row/column labels are emphasised.
 */
export function CrossPcfMatrix({
    result,
    labels,
    r,
    cell,
    labelSize,
    rowLabels,
    colLabels,
    selected,
    focus,
    onFocus,
    onPaint,
}: Props) {
    const view = useMemo(() => matrixAt(result, r), [result, r]);
    const n = view.types.length;
    const ox = rowLabels ? labelSize : 0;
    const oy = colLabels ? labelSize : 0;
    const side = n * cell;
    const label = (t: number) => labels[t] ?? String(t);
    const at = (t: number) => view.types.indexOf(t);
    const rings = selected.map((p) => ({ ...p, i: at(p.a), j: at(p.b) })).filter((p) => p.i >= 0 && p.j >= 0);
    const f = focus ? { i: at(focus.a), j: at(focus.b) } : null;
    const focusCell = f && f.i >= 0 && f.j >= 0 ? f : null;
    const rowColour = new Map<number, string>();
    const colColour = new Map<number, string>();
    for (const p of rings) {
        if (!rowColour.has(p.i)) rowColour.set(p.i, p.colour);
        if (!colColour.has(p.j)) colColour.set(p.j, p.colour);
    }
    const fontSize = Math.min(12, Math.max(cell * 0.8, 8));

    // The stroke lives in a ref: pointer events can arrive faster than renders (a quick drag, a busy
    // page), and each must see the cell the previous one ended on.
    const stroke = useRef<Cell | null>(null);
    const cellAt = (e: PointerEvent<SVGElement>): Cell | null => {
        const box = e.currentTarget.getBoundingClientRect();
        const i = Math.floor((e.clientY - box.top - oy) / cell);
        const j = Math.floor((e.clientX - box.left - ox) / cell);
        return i >= 0 && j >= 0 && i < n && j < n ? { i, j } : null;
    };
    const pairOf = (c: Cell): Pair => ({ a: view.types[c.i], b: view.types[c.j] });

    if (n === 0) return <div className="p-2 opacity-70">No cells of any type inside the domain.</div>;
    const cellBox = (i: number, j: number, inset: number) => ({
        x: j * cell + inset,
        y: i * cell + inset,
        width: Math.max(cell - 0.5 - 2 * inset, 0.5),
        height: Math.max(cell - 0.5 - 2 * inset, 0.5),
    });

    return (
        <svg
            width={ox + side}
            height={oy + side}
            role="img"
            aria-label="cross-PCF matrix"
            style={{ cursor: onPaint ? "pointer" : undefined, userSelect: "none", touchAction: "none" }}
            onPointerDown={(e) => {
                const c = cellAt(e);
                if (!c || !onPaint || e.button !== 0) return;
                try {
                    e.currentTarget.setPointerCapture(e.pointerId); // keep the stroke when it leaves the svg
                } catch {
                    // not an active pointer (a synthetic event): the stroke still works inside the svg
                }
                e.preventDefault();
                stroke.current = c;
                onPaint({ kind: "start", pair: pairOf(c), modifier: e.shiftKey || e.metaKey || e.ctrlKey });
            }}
            onPointerMove={(e) => {
                const c = cellAt(e);
                const from = stroke.current;
                if (from && c && (c.i !== from.i || c.j !== from.j)) {
                    for (const step of cellsAlong(from, c)) onPaint?.({ kind: "cell", pair: pairOf(step) });
                    stroke.current = c;
                }
                const p = c ? pairOf(c) : null;
                if (p?.a !== focus?.a || p?.b !== focus?.b) onFocus?.(p);
            }}
            onPointerUp={() => {
                stroke.current = null;
            }}
            onPointerCancel={() => {
                stroke.current = null;
            }}
            onPointerLeave={() => {
                if (!stroke.current) onFocus?.(null);
            }}
        >
            <g transform={`translate(${ox},${oy})`}>
                {view.types.map((ti, i) =>
                    view.types.map((tj, j) => (
                        <rect
                            key={`${ti}-${tj}`}
                            {...cellBox(i, j, 0)}
                            fill={view.pairs[i * n + j] > 0 ? gColour(view.value[i * n + j]) : "rgb(60,60,60)"}
                        />
                    )),
                )}
                <g pointerEvents="none" fill="none">
                    {rings.map((p) => (
                        <g key={`ring-${p.a}-${p.b}`}>
                            <rect {...cellBox(p.i, p.j, 1)} stroke="black" strokeWidth={3} />
                            <rect {...cellBox(p.i, p.j, 1)} stroke={p.colour} strokeWidth={1.75} />
                        </g>
                    ))}
                    {focusCell && (
                        <g>
                            <rect {...cellBox(focusCell.i, focusCell.j, -2)} stroke="black" strokeWidth={4} />
                            <rect {...cellBox(focusCell.i, focusCell.j, -2)} stroke="white" strokeWidth={2} />
                        </g>
                    )}
                </g>
                {rowLabels &&
                    view.types.map((t, i) => (
                        <text
                            key={`row-${t}`}
                            x={-4}
                            y={i * cell + cell / 2}
                            textAnchor="end"
                            dominantBaseline="middle"
                            fontSize={fontSize}
                            fontWeight={focusCell?.i === i ? 800 : rowColour.has(i) ? 600 : 400}
                            fill={rowColour.get(i) ?? "currentColor"}
                            textDecoration={focusCell?.i === i ? "underline" : undefined}
                        >
                            {label(t)}
                        </text>
                    ))}
                {colLabels &&
                    view.types.map((t, j) => (
                        <text
                            key={`col-${t}`}
                            transform={`translate(${j * cell + cell / 2},-4) rotate(-60)`}
                            fontSize={fontSize}
                            fontWeight={focusCell?.j === j ? 800 : colColour.has(j) ? 600 : 400}
                            fill={colColour.get(j) ?? "currentColor"}
                            textDecoration={focusCell?.j === j ? "underline" : undefined}
                        >
                            {label(t)}
                        </text>
                    ))}
            </g>
        </svg>
    );
}
