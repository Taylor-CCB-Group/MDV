import type { AxisConfig } from "@/react/scatter_state";
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

type DotPlotPanelProps = {
    panel: number;
    title?: string;
    /** Groups for every row category code in this panel. */
    groups: DotPlotGroup[];
    /** Row category codes to draw, top to bottom. */
    rows: number[];
    fieldLabels: string[];
    layout: DotPlotLayout;
    xAxis: AxisConfig;
    fracMax: number;
    meanColor: (mean: number) => string;
    hovered: DotPlotHover | null;
    onHover: (hover: DotPlotHover | null) => void;
};

/** One dot plot: a dot per (row, field), its x axis, and an optional title above. */
export default function DotPlotPanel({
    panel,
    title,
    groups,
    rows,
    fieldLabels,
    layout,
    xAxis,
    fracMax,
    meanColor,
    hovered,
    onHover,
}: DotPlotPanelProps) {
    const left = layout.panelLefts[panel];
    const { top, panelWidth, plotHeight, cellWidth, cellHeight, maxRadius } = layout;
    const xScale = useMemo(
        () => scaleBand({ domain: fieldLabels.map((_, f) => f), range: [0, panelWidth] }),
        [fieldLabels, panelWidth],
    );
    const { tickfont, rotate_labels: rotate } = xAxis;

    return (
        <g transform={`translate(${left},${top})`}>
            {title !== undefined && (
                <text
                    x={panelWidth / 2}
                    y={-6}
                    textAnchor="middle"
                    fontSize={tickfont + 2}
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
                tickFormat={(f) => fieldLabels[f]}
                stroke="var(--text_color)"
                tickStroke="var(--text_color)"
                tickLabelProps={{
                    fill: "var(--text_color)",
                    fontSize: tickfont,
                    textAnchor: rotate ? "end" : "middle",
                    verticalAnchor: rotate ? "middle" : "start",
                    angle: rotate ? -45 : 0,
                }}
            />
        </g>
    );
}
