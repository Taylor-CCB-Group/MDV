import { AxisBottom, AxisLeft } from "@visx/axis";
import { scaleBand, scaleLinear } from "@visx/scale";
import { type PointerEvent, useMemo } from "react";
import type { MultiLineAxis } from "./multiLineConfig";
import type { DensityCurve } from "./multiLineStats";

/** Space kept free above and to the right of the plot, for the last tick labels. */
const PAD_TOP = 10;
const PAD_RIGHT = 15;
/** Short ticks, so that the tick labels and the axis label fit in the legacy chart's 30px. */
const TICK_LENGTH = 4;
const LABEL_OFFSET = 2;
/** Room each x tick label needs. */
const X_TICK_SPACING = 50;

const linePath = (points: [number, number][]) =>
    points.length ? `M${points.map(([x, y]) => `${x},${y}`).join("L")}` : "";

/** The line closed down to `bottom`, for the fill under it. */
const areaPath = (points: [number, number][], bottom: number) =>
    points.length ? `${linePath(points)}L${points[points.length - 1][0]},${bottom}L${points[0][0]},${bottom}Z` : "";

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
    xLabel,
    yLabel,
    onCurveHover,
    onCurveLeave,
}: MultiLinePanelProps) {
    const left = axis.y.size;
    const plotWidth = Math.max(0, width - left - PAD_RIGHT);
    const plotHeight = Math.max(0, height - axis.x.size - PAD_TOP);

    const xScale = useMemo(() => scaleLinear({ domain, range: [0, plotWidth] }), [domain, plotWidth]);
    const rowScale = useMemo(
        () =>
            scaleBand({
                domain: curves.map((c) => c.category),
                range: [0, plotHeight],
            }),
        [curves, plotHeight],
    );

    const paths = useMemo(() => {
        const xs = Array.from(ticks, (t) => xScale(t));
        const rowHeight = stacked ? rowScale.bandwidth() : plotHeight;
        return curves.map((curve) => {
            const bottom = stacked ? (rowScale(curve.category) ?? 0) + rowHeight : plotHeight;
            const scale = rowHeight / curve.max;
            const ys = Array.from(curve.density, (d) => bottom - d * scale);
            const points = ys.map((y, i) => [xs[i], y] as [number, number]);
            return {
                category: curve.category,
                line: linePath(points),
                area: fill ? areaPath(points, bottom) : null,
            };
        });
    }, [curves, ticks, xScale, rowScale, stacked, fill, plotHeight]);

    if (plotWidth === 0 || plotHeight === 0) return null;
    const tickLabelProps = (fontSize: number) => ({ fill: "var(--text_color)", fontSize });
    return (
        <svg width={width} height={height} className="block">
            <g transform={`translate(${left},${PAD_TOP})`}>
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
                <AxisBottom
                    top={plotHeight}
                    scale={xScale}
                    numTicks={Math.max(3, Math.floor(plotWidth / X_TICK_SPACING))}
                    tickLength={TICK_LENGTH}
                    stroke="var(--text_color)"
                    tickStroke="var(--text_color)"
                    tickLabelProps={tickLabelProps(axis.x.tickfont)}
                    label={xLabel}
                    labelOffset={LABEL_OFFSET}
                    labelProps={{ ...tickLabelProps(axis.x.tickfont), textAnchor: "middle" }}
                />
                {stacked ? (
                    <AxisLeft
                        scale={rowScale}
                        tickFormat={(code) => categoryLabels[code] ?? ""}
                        stroke="var(--text_color)"
                        tickStroke="var(--text_color)"
                        tickLabelProps={{ ...tickLabelProps(axis.y.tickfont), textAnchor: "end", dy: "0.33em" }}
                    />
                ) : (
                    // each curve is scaled to its own peak, so a density scale would mean nothing
                    <AxisLeft
                        scale={scaleLinear({ domain: [0, 1], range: [plotHeight, 0] })}
                        numTicks={0}
                        hideTicks
                        stroke="var(--text_color)"
                        label={yLabel}
                        labelOffset={LABEL_OFFSET}
                        labelProps={{ ...tickLabelProps(axis.y.tickfont), textAnchor: "middle" }}
                    />
                )}
            </g>
        </svg>
    );
}
