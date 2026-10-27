import type { XLabelLayout } from "@/charts/dotplot/xAxisLabels";
import { AxisBottom } from "@visx/axis";
import { scaleBand } from "@visx/scale";
import { useMemo } from "react";
import type { DotPlotLayout } from "./dotPlotLayout";
import { getDotRadius } from "./dotPlotLayout";
import type { DotPlotGroup } from "./dotPlotStats";

/** One dot: which panel, row category code and field (index into the fields). */
export type DotPlotDot = { panel: number; row: number; field: number };

/** Tick length of d3's axisBottom, which `layoutXLabels`' spacing assumes. */
const TICK_LENGTH = 6;
/** Outline of the dot whose cells are filtered, as in the legacy dot plot. */
const SELECTED_STROKE_WIDTH = 4;

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
    /** The dot whose cells are filtered, if any. */
    selected: DotPlotDot | null;
    onDotPointer: (event: MouseEvent, dot: DotPlotDot) => void;
    onDotClick: (dot: DotPlotDot) => void;
    onLabelPointer: (event: MouseEvent, field: number) => void;
    onPointerLeave: () => void;
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
    selected,
    onDotPointer,
    onDotClick,
    onLabelPointer,
    onPointerLeave,
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
                    const dot = { panel, row, field };
                    const isSelected = selected?.panel === panel && selected.row === row && selected.field === field;
                    return (
                        <circle
                            // biome-ignore lint/suspicious/noArrayIndexKey: fields have a fixed order
                            key={`${row}-${field}`}
                            cx={field * cellWidth + cellWidth / 2}
                            cy={cy}
                            r={getDotRadius(frac, fracMax, maxRadius)}
                            fill={meanColor(mean)}
                            stroke="black"
                            strokeWidth={isSelected ? SELECTED_STROKE_WIDTH : 1}
                            cursor="pointer"
                            onPointerMove={(e) => onDotPointer(e.nativeEvent, dot)}
                            onPointerLeave={onPointerLeave}
                            onClick={() => onDotClick(dot)}
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
                    onPointerMove: (e) => onLabelPointer(e.nativeEvent, field),
                    onPointerLeave,
                })}
            />
        </g>
    );
}
