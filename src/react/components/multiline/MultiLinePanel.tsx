import { scaleBand, scaleLinear } from "@visx/scale";
import { type PointerEvent, useMemo } from "react";
import { canvasMeasureText } from "@/charts/dotplot/xAxisLabels";
import type { MultiLineAxis } from "./multiLineConfig";
import { AXIS_FONT_FAMILY, getMultiLineLayout, LABEL_GAP, LINE_HEIGHT, TICK_LENGTH } from "./multiLineLayout";
import type { DensityCurve } from "./multiLineStats";

const linePath = (points: [number, number][]) =>
    points.length ? `M${points.map(([x, y]) => `${x},${y}`).join("L")}` : "";

/** The line closed down to `bottom`, for the fill under it. */
const areaPath = (points: [number, number][], bottom: number) =>
    points.length ? `${linePath(points)}L${points[points.length - 1][0]},${bottom}L${points[0][0]},${bottom}Z` : "";

const TEXT = "var(--text_color)";

export type MultiLinePanelProps = {
    width: number;
    height: number;
    /** x domain, in data units. */
    domain: [number, number];
    /** Where each curve's density was evaluated, in data units. */
    ticks: ArrayLike<number>;
    curves: DensityCurve[];
    /** Label of each category code. */
    categoryLabels: readonly string[];
    /** Colour of each category code. */
    colors: readonly string[];
    stacked: boolean;
    fill: boolean;
    axis: { x: MultiLineAxis; y: MultiLineAxis };
    /** Fit the margins to the labels; off, the axis sizes are used. */
    fitLabels: boolean;
    xLabel: string;
    /** Shown beside the y axis when not stacked; stacked, the rows are labelled instead. */
    yLabel: string;
    onCurveHover?: (event: PointerEvent<SVGPathElement>, category: number) => void;
    onCurveLeave?: () => void;
};

/**
 * Density curves for one value column, one per category. Each curve is scaled to its own peak,
 * either over the whole plot or, stacked, within its own row. Draws from data alone: no loading,
 * settings or chart, so other charts can hold several.
 */
export default function MultiLinePanel({
    width,
    height,
    domain,
    ticks,
    curves,
    categoryLabels,
    colors,
    stacked,
    fill,
    axis,
    fitLabels,
    xLabel,
    yLabel,
    onCurveHover,
    onCurveLeave,
}: MultiLinePanelProps) {
    const fullRowLabels = useMemo(
        () => (stacked ? curves.map((c) => categoryLabels[c.category] ?? "") : null),
        [stacked, curves, categoryLabels],
    );
    const { size: xSize, tickfont: xFont } = axis.x;
    const { size: ySize, tickfont: yFont } = axis.y;
    const layout = useMemo(
        () =>
            getMultiLineLayout(
                {
                    width,
                    height,
                    // axis is mutated in place by its settings, so it's rebuilt from its values
                    axis: { x: { size: xSize, tickfont: xFont }, y: { size: ySize, tickfont: yFont } },
                    fitLabels,
                    domain,
                    xLabel,
                    yLabel,
                    rowLabels: fullRowLabels,
                },
                canvasMeasureText,
            ),
        [width, height, xSize, xFont, ySize, yFont, fitLabels, domain, xLabel, yLabel, fullRowLabels],
    );
    const { left, top, plotWidth, plotHeight } = layout;

    const xScale = useMemo(() => scaleLinear({ domain, range: [0, plotWidth] }), [domain, plotWidth]);
    const rowScale = useMemo(
        () => scaleBand({ domain: curves.map((c) => c.category), range: [0, plotHeight] }),
        [curves, plotHeight],
    );

    const paths = useMemo(() => {
        const xs = Array.from(ticks, (t) => xScale(t));
        const rowHeight = stacked ? rowScale.bandwidth() : plotHeight;
        return curves.map((curve) => {
            const bottom = stacked ? (rowScale(curve.category) ?? 0) + rowHeight : plotHeight;
            const scale = rowHeight / curve.max;
            const points = Array.from(curve.density, (d, i) => [xs[i], bottom - d * scale] as [number, number]);
            return {
                category: curve.category,
                line: linePath(points),
                area: fill ? areaPath(points, bottom) : null,
            };
        });
    }, [curves, ticks, xScale, rowScale, stacked, fill, plotHeight]);

    if (plotWidth === 0 || plotHeight === 0) return null;
    const tickLabelTop = TICK_LENGTH + LABEL_GAP;
    const { rowLabels } = layout;
    return (
        <svg width={width} height={height} className="block" fontFamily={AXIS_FONT_FAMILY}>
            <g transform={`translate(${left},${top})`}>
                {paths.map(
                    (p) =>
                        p.area && (
                            <path
                                key={`area-${p.category}`}
                                d={p.area}
                                fill={colors[p.category]}
                                opacity={0.5}
                                pointerEvents="none"
                            />
                        ),
                )}
                {paths.map((p) => (
                    <path
                        key={p.category}
                        d={p.line}
                        fill="none"
                        stroke={colors[p.category]}
                        strokeWidth={2}
                        strokeLinejoin="round"
                        onPointerMove={onCurveHover && ((e) => onCurveHover(e, p.category))}
                        onPointerLeave={onCurveLeave}
                    />
                ))}

                <g transform={`translate(0,${plotHeight})`} fill={TEXT} fontSize={xFont}>
                    <line x2={plotWidth} stroke={TEXT} />
                    {layout.xTicks.map((t) => (
                        <g key={t} transform={`translate(${xScale(t)},0)`}>
                            <line y2={TICK_LENGTH} stroke={TEXT} />
                            <text y={tickLabelTop} textAnchor="middle" dominantBaseline="hanging">
                                {layout.xTickFormat(t)}
                            </text>
                        </g>
                    ))}
                    <text
                        x={plotWidth / 2}
                        y={tickLabelTop + xFont * LINE_HEIGHT + LABEL_GAP}
                        textAnchor="middle"
                        dominantBaseline="hanging"
                    >
                        {layout.xLabel !== xLabel && <title>{xLabel}</title>}
                        {layout.xLabel}
                    </text>
                </g>

                <g fill={TEXT}>
                    <line y2={plotHeight} stroke={TEXT} />
                    {rowLabels && fullRowLabels ? (
                        curves.map((c, i) => (
                            <g
                                key={c.category}
                                transform={`translate(0,${(rowScale(c.category) ?? 0) + rowScale.bandwidth() / 2})`}
                            >
                                <line x2={-TICK_LENGTH} stroke={TEXT} />
                                <text
                                    x={-TICK_LENGTH - LABEL_GAP}
                                    textAnchor="end"
                                    dominantBaseline="central"
                                    fontSize={layout.rowFontSize}
                                >
                                    {rowLabels[i] !== fullRowLabels[i] && <title>{fullRowLabels[i]}</title>}
                                    {rowLabels[i]}
                                </text>
                            </g>
                        ))
                    ) : (
                        // each curve is scaled to its own peak, so a density scale would mean nothing
                        <text
                            transform={`translate(${-LABEL_GAP - 1},${plotHeight / 2}) rotate(-90)`}
                            textAnchor="middle"
                            fontSize={yFont}
                        >
                            {layout.yLabel}
                        </text>
                    )}
                </g>
            </g>
        </svg>
    );
}
