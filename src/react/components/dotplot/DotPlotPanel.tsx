import type { XLabelLayout } from "@/charts/dotplot/xAxisLabels";
import { AxisBottom } from "@visx/axis";
import { scaleBand } from "@visx/scale";
import { useMemo } from "react";
import type { DotPlotLayout } from "./dotPlotLayout";
import { getDotRadius } from "./dotPlotLayout";
import type { DotPlotGroup } from "./dotPlotStats";

export type DotPlotHover = {
    panel: number;
    /** Row category code. */
    row: number;
    /** Index into the fields. */
    field: number;
    /** Dot centre, relative to the chart. */
    x: number;
    y: number;
};

export type DotPlotLabelHover = {
    /** Index into the fields. */
    field: number;
    /** Label's tick, relative to the chart. */
    x: number;
    y: number;
};

/** Tick length of d3's axisBottom, which `layoutXLabels`' spacing assumes. */
const TICK_LENGTH = 6;

type DotPlotPanelProps = {
    panel: number;
    title?: string;
    /** Groups for every row category code in this panel. */
    groups: DotPlotGroup[];
    /** Row category codes to draw, top to bottom. */
    rows: number[];
    xLabels: XLabelLayout;
    layout: DotPlotLayout;
    /** Font size for the panel title. */
    titleFontSize: number;
    fracMax: number;
    meanColor: (mean: number) => string;
    hovered: DotPlotHover | null;
    onHover: (hover: DotPlotHover | null) => void;
    onLabelHover: (hover: DotPlotLabelHover | null) => void;
};

/** One dot plot: a dot per (row, field), its x axis, and an optional title above. */
export default function DotPlotPanel({
    panel,
    title,
    groups,
    rows,
    xLabels,
    layout,
    titleFontSize,
    fracMax,
    meanColor,
    hovered,
    onHover,
    onLabelHover,
}: DotPlotPanelProps) {
    const left = layout.panelLefts[panel];
    const { top, panelWidth, plotHeight, cellWidth, cellHeight, maxRadius } = layout;
    const rotated = xLabels.angle === 45;
    const fieldCount = xLabels.labels.length;
    const xScale = useMemo(
        () => scaleBand({ domain: Array.from({ length: fieldCount }, (_, f) => f), range: [0, panelWidth] }),
        [fieldCount, panelWidth],
    );

    return (
        <g transform={`translate(${left},${top})`}>
            {title !== undefined && (
                <text
                    x={panelWidth / 2}
                    y={-6}
                    textAnchor="middle"
                    fontSize={titleFontSize}
                    fontWeight={500}
                    fill="var(--text_color)"
                >
                    {title}
                </text>
            )}
            <rect width={panelWidth} height={plotHeight} fill="none" stroke="var(--text_color)" strokeOpacity={0.2} />
            {rows.map((row, i) => {
                const cells = groups[row]?.cells ?? [];
                const cy = i * cellHeight + cellHeight / 2;
                return cells.map(({ mean, frac }, field) => {
                    const cx = field * cellWidth + cellWidth / 2;
                    const isHovered = hovered?.panel === panel && hovered.row === row && hovered.field === field;
                    return (
                        <circle
                            // biome-ignore lint/suspicious/noArrayIndexKey: fields have a fixed order
                            key={`${row}-${field}`}
                            cx={cx}
                            cy={cy}
                            r={getDotRadius(frac, fracMax, maxRadius)}
                            fill={meanColor(mean)}
                            stroke="black"
                            strokeWidth={isHovered ? 2 : 1}
                            onPointerEnter={() => onHover({ panel, row, field, x: left + cx, y: top + cy })}
                            onPointerLeave={() => onHover(null)}
                        />
                    );
                });
            })}
            <AxisBottom
                top={plotHeight}
                scale={xScale}
                tickLength={TICK_LENGTH}
                tickFormat={(field) => xLabels.labels[field]}
                stroke="var(--text_color)"
                tickStroke="var(--text_color)"
                tickLabelProps={(field) => ({
                    fill: "var(--text_color)",
                    fontSize: xLabels.fontSize,
                    textAnchor: rotated ? "end" : "middle",
                    verticalAnchor: rotated ? "middle" : "start",
                    angle: rotated ? -45 : 0,
                    onPointerEnter: () =>
                        onLabelHover({
                            field,
                            x: left + field * cellWidth + cellWidth / 2,
                            y: top + plotHeight + TICK_LENGTH,
                        }),
                    onPointerLeave: () => onLabelHover(null),
                })}
            />
        </g>
    );
}
