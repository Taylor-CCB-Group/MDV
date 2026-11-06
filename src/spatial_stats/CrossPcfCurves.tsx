import type { PcfMatrixBinnedResult } from "intraspatial/spatial/pcf";
import { useId } from "react";
import { pairCurve, type RSelection } from "./CrossPcfMatrix";
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
    r: RSelection;
    unit: string;
    width: number;
    height: number;
};

const niceMax = (v: number) => (v <= 2 ? 2 : v <= 3 ? 3 : v <= 5 ? 5 : Math.ceil(v / 5) * 5);

/** g(r) for the chosen pairs, with the g = 1 line and the matrix's current r bin shaded. */
export function CrossPcfCurves({ result, pairs, preview, focus, onFocus, r, unit, width, height }: Props) {
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
    const path = (g: number[]) =>
        g
            .map((v, k) => `${k ? "L" : "M"}${sx(centres[k]).toFixed(1)},${sy(Number.isFinite(v) ? v : 0).toFixed(1)}`)
            .join("");

    return (
        <svg width={width} height={height} role="img" aria-label="cross-PCF curves" fontSize={10}>
            <clipPath id={clip}>
                <rect x={left} y={top} width={w} height={h} />
            </clipPath>
            {r.kind === "bin" && (
                <rect
                    x={sx(r.k * dr)}
                    y={top}
                    width={sx(dr) - sx(0)}
                    height={h}
                    fill="rgb(255,200,0)"
                    fillOpacity={0.12}
                />
            )}
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
                    const faded = focus !== null && !focused;
                    return (
                        <g key={c.key}>
                            <path
                                d={path(c.g)}
                                fill="none"
                                stroke={c.colour}
                                strokeWidth={focused ? 3.5 : 1.75}
                                strokeOpacity={faded ? 0.15 : 0.95}
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
    pairs,
    preview,
    focus,
    onFocus,
    onRemove,
    height,
}: {
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
                            opacity: focus && !focused ? 0.45 : 1,
                            fontWeight: focused ? 700 : 400,
                        }}
                        onPointerEnter={() => onFocus?.({ a: p.a, b: p.b })}
                        onPointerLeave={() => onFocus?.(null)}
                    >
                        <span
                            aria-hidden
                            style={{
                                width: 14,
                                height: 0,
                                borderTop: `${focused ? 3 : 2}px ${isPreview ? "dashed" : "solid"} ${p.colour}`,
                                flex: "none",
                            }}
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
