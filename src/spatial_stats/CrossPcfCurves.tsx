import type { PcfMatrixBinnedResult } from "intraspatial/spatial/pcf";
import { type PointerEvent, useId, useRef, useState } from "react";
import { pairCurve } from "./CrossPcfMatrix";
import { nearestBin, type RBin } from "./rBins";
import { useTweened } from "./useTweened";

/** One curve per chosen pair; colours cycle. */
export const PAIR_COLOURS = ["#e8743b", "#4f8fe0", "#2bb3a6", "#c95fd1", "#d4b020", "#e0526f", "#7dbb42", "#8d7ae6"];

export type CurvePair = { a: number; b: number; label: string; colour: string };

type Pair = { a: number; b: number };
const isPair = (p: Pair, q: Pair | null) => !!q && p.a === q.a && p.b === q.b;

type Props = {
    result: PcfMatrixBinnedResult;
    pairs: CurvePair[];
    /** A transient extra curve (a focused pair that isn't chosen), drawn dashed. */
    preview?: CurvePair;
    /** The pair highlighted across the chart: its curve is drawn on top and the others fade. */
    focus: Pair | null;
    onFocus?: (pair: Pair | null) => void;
    /** The r's shown, each with its band colour; the mean has no band. */
    rs: readonly { bin: RBin; colour: string }[];
    /** Dragging in the plot moves the nearest band (r `index`) to the bin under the pointer. */
    onMoveR?: (index: number, k: number) => void;
    /** Shift-press in the plot (or any press when only the mean is shown) adds an r there. */
    onAddR?: (k: number) => void;
    /** The r (index into `rs`) highlighted across the chart; hovering or dragging a band sets it. */
    rFocus: number | null;
    onRFocus?: (index: number | null) => void;
    unit: string;
    width: number;
    height: number;
};

/** The r bin under plot x-coordinate `x` (plot spans `[left, left + w]` for r ∈ [0, rSpan)), clamped. */
export function binAtX(x: number, left: number, w: number, rSpan: number, dr: number, nBins: number): number {
    const rv = ((x - left) / Math.max(w, 1e-9)) * rSpan;
    return Math.min(nBins - 1, Math.max(0, Math.floor(rv / dr)));
}

/** A g(r) curve as an SVG path in a `width × height` box over `[0, yMax]`; values above `yMax` are
 *  clipped to the top. */
export function miniCurvePath(g: readonly number[], width: number, height: number, yMax: number): string {
    const n = g.length;
    return g
        .map((v, k) => {
            const x = n > 1 ? (k / (n - 1)) * width : 0;
            const y = height - (Math.min(Number.isFinite(v) ? v : 0, yMax) / yMax) * height;
            return `${k ? "L" : "M"}${x.toFixed(1)},${y.toFixed(1)}`;
        })
        .join("");
}

const niceMax = (v: number) => (v <= 2 ? 2 : v <= 3 ? 3 : v <= 5 ? 5 : Math.ceil(v / 5) * 5);

/** g(r) for the chosen pairs, with the g = 1 line and the matrix's current r bin shaded. */
export function CrossPcfCurves({
    result,
    pairs,
    preview,
    focus,
    onFocus,
    rs,
    onMoveR,
    onAddR,
    rFocus,
    onRFocus,
    unit,
    width,
    height,
}: Props) {
    const left = 34;
    const right = 8;
    const top = 8;
    const bottom = 30;
    const w = Math.max(0, width - left - right);
    const h = Math.max(0, height - top - bottom);
    const { r: centres, nBins: B } = result;
    const dr = centres[0] * 2;
    const rMax = dr * B;
    const shown = preview ? [...pairs, preview] : pairs;
    const curves = shown.map((p) => ({
        ...p,
        key: `${p.a}-${p.b}${p === preview ? "-preview" : ""}`,
        g: pairCurve(result, p.a, p.b),
    }));
    let gMax = 2;
    for (const c of curves) for (const v of c.g) if (Number.isFinite(v)) gMax = Math.max(gMax, v);
    const yMax = niceMax(gMax);
    // Axis ranges ease to their new values; ticks are the new scale's, sliding into place.
    const yShown = useTweened(yMax);
    const rShown = useTweened(rMax);
    const clip = useId();
    const sx = (rv: number) => left + (rv / rShown) * w;
    const sy = (gv: number) => top + h - (gv / yShown) * h;
    const inY = (t: number) => t <= yShown * 1.0001;
    const inX = (t: number) => t <= rShown * 1.0001;
    const yTicks = Array.from({ length: Math.min(yMax, 5) + 1 }, (_, i) => (i * yMax) / Math.min(yMax, 5)).filter(inY);
    const xTicks = Array.from({ length: 6 }, (_, i) => (i * rMax) / 5).filter(inX);
    // The band being dragged (and the bin it is on) lives in a ref, so moves arriving before the next
    // render see the stroke; the state copy only drives drawing.
    const drag = useRef<{ index: number; k: number } | null>(null);
    const [dragIndex, setDragIndex] = useState<number | null>(null);
    const bins = rs.map((x) => x.bin);
    const endDrag = () => {
        drag.current = null;
        setDragIndex(null);
    };
    const binAt = (e: PointerEvent<SVGSVGElement>) => {
        const box = e.currentTarget.getBoundingClientRect();
        return binAtX(e.clientX - box.left, left, w, rShown, dr, B);
    };
    const inPlot = (e: PointerEvent<SVGSVGElement>) => {
        const box = e.currentTarget.getBoundingClientRect();
        const x = e.clientX - box.left;
        const y = e.clientY - box.top;
        return x >= left && x <= left + w && y >= top && y <= top + h + bottom;
    };
    const path = (g: number[]) =>
        g
            .map((v, k) => `${k ? "L" : "M"}${sx(centres[k]).toFixed(1)},${sy(Number.isFinite(v) ? v : 0).toFixed(1)}`)
            .join("");

    return (
        <svg
            width={width}
            height={height}
            role="img"
            aria-label="cross-PCF curves"
            fontSize={10}
            style={{
                cursor: dragIndex !== null ? "ew-resize" : onMoveR ? "col-resize" : undefined,
                userSelect: "none",
                touchAction: "none",
            }}
            onPointerDown={(e) => {
                if (!onMoveR || !onAddR || e.button !== 0 || !inPlot(e)) return;
                try {
                    e.currentTarget.setPointerCapture(e.pointerId);
                } catch {
                    // a synthetic pointer: the drag still works inside the svg
                }
                e.preventDefault();
                const k = binAt(e);
                let index = nearestBin(bins, k);
                if (e.shiftKey || index < 0) {
                    // add a band here (or take hold of the one already on this bin) and keep dragging it
                    const existing = bins.indexOf(k);
                    index = existing >= 0 ? existing : bins.length;
                    if (existing < 0) onAddR(k);
                } else if (bins[index] !== k) onMoveR(index, k);
                drag.current = { index, k };
                setDragIndex(index);
                onRFocus?.(index);
            }}
            onPointerMove={(e) => {
                const d = drag.current;
                if (!d) {
                    // not dragging: the band under the pointer, if any, is the focused r
                    const over = inPlot(e) ? bins.indexOf(binAt(e)) : -1;
                    const next = over >= 0 ? over : null;
                    if (next !== rFocus) onRFocus?.(next);
                    return;
                }
                if (!onMoveR) return;
                const k = binAt(e);
                if (k !== d.k) {
                    onMoveR(d.index, k);
                    drag.current = { ...d, k };
                }
            }}
            onPointerUp={endDrag}
            onPointerCancel={endDrag}
            onPointerLeave={() => {
                if (!drag.current && rFocus !== null) onRFocus?.(null);
            }}
        >
            <clipPath id={clip}>
                <rect x={left} y={top} width={w} height={h} />
            </clipPath>
            {rs.map(({ bin, colour }, index) => {
                if (bin === "mean") return null;
                const focused = rFocus === index || dragIndex === index;
                const receded = rFocus !== null && !focused;
                const x = sx(bin * dr);
                const bw = sx(dr) - sx(0);
                return (
                    <g key={`band-${bin}`} pointerEvents="none">
                        <rect
                            x={x}
                            y={top}
                            width={bw}
                            height={h}
                            fill={colour}
                            fillOpacity={focused ? 0.3 : receded ? 0.06 : 0.13}
                        />
                        {focused && (
                            <>
                                <line x1={x} x2={x} y1={top} y2={top + h} stroke={colour} strokeWidth={1.5} />
                                <line x1={x + bw} x2={x + bw} y1={top} y2={top + h} stroke={colour} strokeWidth={1.5} />
                            </>
                        )}
                        {/* grab handle on the axis, so the band reads as draggable */}
                        <path
                            d={`M${sx((bin + 0.5) * dr)},${top + h + 1} ${focused ? "l-7,10 h14" : "l-5,7 h10"} z`}
                            fill={colour}
                            fillOpacity={receded ? 0.5 : 1}
                            stroke="black"
                            strokeWidth={0.75}
                        />
                    </g>
                );
            })}
            {yTicks.map((t) => (
                <g key={`y${t}`}>
                    <line x1={left} x2={left + w} y1={sy(t)} y2={sy(t)} stroke="currentColor" strokeOpacity={0.12} />
                    <text
                        x={left - 4}
                        y={sy(t)}
                        textAnchor="end"
                        dominantBaseline="middle"
                        fill="currentColor"
                        opacity={0.7}
                    >
                        {t}
                    </text>
                </g>
            ))}
            <line
                x1={left}
                x2={left + w}
                y1={sy(1)}
                y2={sy(1)}
                stroke="currentColor"
                strokeOpacity={0.5}
                strokeDasharray="4 3"
            />
            {xTicks.map((t) => (
                <text key={`x${t}`} x={sx(t)} y={top + h + 12} textAnchor="middle" fill="currentColor" opacity={0.7}>
                    {Math.round(t)}
                </text>
            ))}
            <text x={left + w / 2} y={top + h + 26} textAnchor="middle" fill="currentColor" opacity={0.8}>
                {`r (${unit})`}
            </text>
            <g clipPath={`url(#${clip})`}>
                {/* Focused curve last, so it is drawn on top. */}
                {[...curves.filter((c) => !isPair(c, focus)), ...curves.filter((c) => isPair(c, focus))].map((c) => {
                    const focused = isPair(c, focus);
                    // Hovering an unchosen pair (the dashed preview) only recedes the chosen curves a little,
                    // so they stay comparable with it; focusing a chosen one dims the rest more.
                    const fade = focus === null || focused ? 0.95 : preview ? 0.75 : 0.4;
                    return (
                        <g key={c.key}>
                            <path
                                d={path(c.g)}
                                fill="none"
                                stroke={c.colour}
                                strokeWidth={focused ? 3.5 : 1.75}
                                strokeOpacity={fade}
                                strokeDasharray={c.key.endsWith("-preview") ? "5 3" : undefined}
                                pointerEvents="none"
                            />
                            {/* a wide invisible line makes thin curves easy to hover */}
                            <path
                                d={path(c.g)}
                                fill="none"
                                stroke="transparent"
                                strokeWidth={10}
                                pointerEvents="stroke"
                                style={{ cursor: "pointer" }}
                                onPointerEnter={() => onFocus?.({ a: c.a, b: c.b })}
                                onPointerLeave={() => onFocus?.(null)}
                            />
                        </g>
                    );
                })}
            </g>
            {curves.length === 0 && (
                <text x={left + w / 2} y={top + h / 2} textAnchor="middle" fill="currentColor" opacity={0.6}>
                    Drag across matrix cells to choose A → B pairs
                </text>
            )}
        </svg>
    );
}

/** The chosen pairs as a list beside the plot: hover highlights the pair everywhere, × removes it. */
export function PairLegend({
    result,
    pairs,
    preview,
    focus,
    onFocus,
    onRemove,
    height,
}: {
    result: PcfMatrixBinnedResult;
    pairs: CurvePair[];
    preview?: CurvePair;
    focus: Pair | null;
    onFocus?: (pair: Pair | null) => void;
    onRemove?: (pair: Pair) => void;
    height: number;
}) {
    const items = preview ? [...pairs, preview] : pairs;
    return (
        <ul className="m-0 p-0 overflow-y-auto text-xs" style={{ listStyle: "none", maxHeight: height, minWidth: 96 }}>
            {items.map((p) => {
                const focused = isPair(p, focus);
                const isPreview = p === preview;
                return (
                    <li
                        key={`${p.a}-${p.b}${isPreview ? "-preview" : ""}`}
                        className="flex items-center gap-1 px-1 rounded"
                        style={{
                            background: focused ? "rgba(255,255,255,0.12)" : undefined,
                            opacity: focus === null || focused ? 1 : preview ? 0.8 : 0.55,
                            fontWeight: focused ? 700 : 400,
                        }}
                        onPointerEnter={() => onFocus?.({ a: p.a, b: p.b })}
                        onPointerLeave={() => onFocus?.(null)}
                    >
                        <MiniCurve
                            g={pairCurve(result, p.a, p.b)}
                            colour={p.colour}
                            dashed={isPreview}
                            bold={focused}
                        />
                        <span className="truncate flex-1" style={{ color: p.colour }}>
                            {p.label}
                        </span>
                        {!isPreview && onRemove && (
                            <button
                                type="button"
                                aria-label={`Remove ${p.label}`}
                                className="opacity-60 hover:opacity-100 px-1"
                                style={{ background: "none", border: "none", color: "inherit", cursor: "pointer" }}
                                onClick={() => onRemove({ a: p.a, b: p.b })}
                            >
                                ×
                            </button>
                        )}
                    </li>
                );
            })}
        </ul>
    );
}

/** The legend's swatch: the pair's own g(r) in miniature, on its own scale (at least [0, 2], so a
 *  flat g ≈ 1 sits mid-box), with a faint g = 1 line. */
function MiniCurve({ g, colour, dashed, bold }: { g: number[]; colour: string; dashed: boolean; bold: boolean }) {
    const w = 30;
    const h = 14;
    let max = 0;
    for (const v of g) if (Number.isFinite(v)) max = Math.max(max, v);
    const yMax = niceMax(max);
    const one = h - (1 / yMax) * h;
    return (
        <svg aria-hidden width={w} height={h} style={{ flex: "none", overflow: "visible" }}>
            <line x1={0} x2={w} y1={one} y2={one} stroke="currentColor" strokeOpacity={0.25} strokeWidth={0.75} />
            <path
                d={miniCurvePath(g, w, h, yMax)}
                fill="none"
                stroke={colour}
                strokeWidth={bold ? 2.25 : 1.5}
                strokeDasharray={dashed ? "3 2" : undefined}
            />
        </svg>
    );
}
