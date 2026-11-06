import type { PcfMatrixBinnedResult } from "intraspatial/spatial/pcf";
import { useMemo, useState } from "react";

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

/** A chosen pair as stored in config: category labels, anchor A → neighbour B (g is asymmetric). */
export type LabelPair = { a: string; b: string };

/** Configured pairs as type ids of this result; pairs naming a category it doesn't have are dropped. */
export function resolvePairs(labels: readonly string[], pairs: readonly LabelPair[]) {
    const out: { a: number; b: number; index: number }[] = [];
    pairs.forEach((p, index) => {
        const a = labels.indexOf(p.a);
        const b = labels.indexOf(p.b);
        if (a >= 0 && b >= 0) out.push({ a, b, index });
    });
    return out;
}

/** `pairs` with A → B added, or removed if already there. */
export function togglePair(pairs: readonly LabelPair[], a: string, b: string): LabelPair[] {
    const has = pairs.some((p) => p.a === a && p.b === b);
    return has ? pairs.filter((p) => !(p.a === a && p.b === b)) : [...pairs, { a, b }];
}

/** Diverging colour on log2 g: blue = avoidance (g < 1), white = independence, red = attraction. */
export function gColour(g: number, maxLog2 = 2): string {
    if (!Number.isFinite(g) || g <= 0) return "rgb(40,40,40)";
    const t = Math.max(-1, Math.min(1, Math.log2(g) / maxLog2));
    const a = Math.abs(t);
    const mix = (c: number) => Math.round(255 + (c - 255) * a);
    return t >= 0 ? `rgb(${mix(200)},${mix(30)},${mix(40)})` : `rgb(${mix(30)},${mix(80)},${mix(200)})`;
}

type Props = {
    result: PcfMatrixBinnedResult;
    labels: string[];
    r: RSelection;
    unit: string;
    width: number;
    height: number;
    /** Chosen pairs as type ids, with their curve colours — owned by the caller (chart config). */
    selected: readonly { a: number; b: number; colour: string }[];
    onToggle?: (a: number, b: number) => void;
    onHover?: (pair: { a: number; b: number } | null) => void;
};

/**
 * K×K heatmap of cross-PCF at one r (row = anchor type A, column = neighbour type B). Clicking a cell
 * adds or removes that pair; chosen pairs are outlined in their curve colour.
 */
export function CrossPcfMatrix({ result, labels, r, unit, width, height, selected, onToggle, onHover }: Props) {
    const view = useMemo(() => matrixAt(result, r), [result, r]);
    const [hover, setHoverState] = useState<{ i: number; j: number } | null>(null);
    const setHover = (h: { i: number; j: number } | null) => {
        setHoverState(h);
        onHover?.(h ? { a: view.types[h.i], b: view.types[h.j] } : null);
    };
    const n = view.types.length;
    const labelW = 70;
    const footH = 18;
    const side = Math.max(0, Math.min(width - labelW - 8, height - labelW - footH));
    const cell = n > 0 ? side / n : 0;
    const label = (t: number) => labels[t] ?? String(t);
    const outline = new Map(selected.map((p) => [`${p.a}-${p.b}`, p.colour]));

    if (n === 0) return <div className="p-2 opacity-70">No cells of any type inside the domain.</div>;
    const rLabel = r.kind === "mean" ? "mean over r" : `r = ${result.r[r.k].toFixed(1)} ${unit}`;

    return (
        <svg width={width} height={height} role="img" aria-label="cross-PCF matrix" onMouseLeave={() => setHover(null)}>
            <g transform={`translate(${labelW},${labelW})`}>
                {view.types.map((ti, i) =>
                    view.types.map((tj, j) => {
                        const ring = outline.get(`${ti}-${tj}`);
                        return (
                            // biome-ignore lint/a11y/noStaticElementInteractions: matrix cells, pointer only for now
                            <rect
                                key={`${ti}-${tj}`}
                                x={j * cell}
                                y={i * cell}
                                width={Math.max(cell - 0.5, 0.5)}
                                height={Math.max(cell - 0.5, 0.5)}
                                fill={view.pairs[i * n + j] > 0 ? gColour(view.value[i * n + j]) : "rgb(60,60,60)"}
                                stroke={ring ?? (hover?.i === i && hover?.j === j ? "currentColor" : "none")}
                                strokeWidth={ring ? 2.5 : 1}
                                style={{ cursor: onToggle ? "pointer" : undefined }}
                                onMouseEnter={() => setHover({ i, j })}
                                onClick={() => onToggle?.(ti, tj)}
                            />
                        );
                    }),
                )}
                {view.types.map((t, i) => (
                    <text
                        key={`row-${t}`}
                        x={-4}
                        y={i * cell + cell / 2}
                        textAnchor="end"
                        dominantBaseline="middle"
                        fontSize={Math.min(11, cell)}
                        fill="currentColor"
                    >
                        {label(t)}
                    </text>
                ))}
                {view.types.map((t, j) => (
                    <text
                        key={`col-${t}`}
                        transform={`translate(${j * cell + cell / 2},-4) rotate(-60)`}
                        fontSize={Math.min(11, cell)}
                        fill="currentColor"
                    >
                        {label(t)}
                    </text>
                ))}
            </g>
            <text
                x={labelW}
                y={labelW + side + 4}
                dominantBaseline="hanging"
                fontSize={11}
                fill="currentColor"
                opacity={0.8}
            >
                {hover
                    ? `${label(view.types[hover.i])} → ${label(view.types[hover.j])}: g = ${view.value[hover.i * n + hover.j].toFixed(2)} (${rLabel}), ${view.pairs[hover.i * n + hover.j]} pairs`
                    : `${rLabel}; row = A (anchor), column = B (neighbour); click to add/remove a pair`}
            </text>
        </svg>
    );
}
