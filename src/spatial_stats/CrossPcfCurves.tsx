import type { PcfMatrixBinnedResult } from "intraspatial/spatial/pcf";
import { pairCurve, type RSelection } from "./CrossPcfMatrix";

/** One curve per chosen pair; colours cycle. */
export const PAIR_COLOURS = ["#e8743b", "#4f8fe0", "#2bb3a6", "#c95fd1", "#d4b020", "#e0526f", "#7dbb42", "#8d7ae6"];

export type CurvePair = { a: number; b: number; label: string; colour: string };

type Props = {
    result: PcfMatrixBinnedResult;
    pairs: CurvePair[];
    /** A transient extra curve (the hovered matrix cell), drawn dashed. */
    preview?: CurvePair;
    r: RSelection;
    unit: string;
    width: number;
    height: number;
};

const niceMax = (v: number) => (v <= 2 ? 2 : v <= 3 ? 3 : v <= 5 ? 5 : Math.ceil(v / 5) * 5);

/** g(r) for the chosen pairs, with the g = 1 line and the matrix's current r bin shaded. */
export function CrossPcfCurves({ result, pairs, preview, r, unit, width, height }: Props) {
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
    const sx = (rv: number) => left + (rv / rMax) * w;
    const sy = (gv: number) => top + h - (Math.min(gv, yMax) / yMax) * h;
    const yTicks = Array.from({ length: Math.min(yMax, 5) + 1 }, (_, i) => (i * yMax) / Math.min(yMax, 5));
    const xTicks = Array.from({ length: 6 }, (_, i) => (i * rMax) / 5);
    const path = (g: number[]) =>
        g
            .map((v, k) => `${k ? "L" : "M"}${sx(centres[k]).toFixed(1)},${sy(Number.isFinite(v) ? v : 0).toFixed(1)}`)
            .join("");

    return (
        <svg width={width} height={height} role="img" aria-label="cross-PCF curves" fontSize={10}>
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
            {curves.map((c) => (
                <path
                    key={c.key}
                    d={path(c.g)}
                    fill="none"
                    stroke={c.colour}
                    strokeWidth={2}
                    strokeDasharray={c === curves.at(-1) && preview ? "5 3" : undefined}
                />
            ))}
            {curves.length === 0 && (
                <text x={left + w / 2} y={top + h / 2} textAnchor="middle" fill="currentColor" opacity={0.6}>
                    Click matrix cells to add A → B pairs
                </text>
            )}
            {/* Legend, top right, one line per curve. */}
            <g transform={`translate(${left + w - 4},${top + 4})`}>
                {curves.map((c, i) => (
                    <text key={`l${c.key}`} y={i * 12} textAnchor="end" dominantBaseline="hanging" fill={c.colour}>
                        {c.label}
                        {c === curves.at(-1) && preview ? " (hover)" : ""}
                    </text>
                ))}
            </g>
        </svg>
    );
}
