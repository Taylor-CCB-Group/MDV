import { useLayoutEffect, useMemo, useRef, useState } from "react";
import { axisBottom } from "d3";
import { scaleLinear } from "d3-scale";
import { select } from "d3-selection";
import { useBrushX } from "@/react/components/histogram/useBrushX";
import {
    DEFAULT_CONTINUOUS_LEGEND_HEIGHT,
    DEFAULT_CONTINUOUS_LEGEND_WIDTH,
} from "@/react/legend/shared/legendConstants";
import type { LegendContinuousSvgProps } from "@/react/legend/shared/legendTypes";
import {
    formatContinuousTick,
    formatLegendLabel,
    getContinuousLegendLayout,
    getGradientStops,
} from "@/react/legend/shared/legendUtils";

export {
    DEFAULT_CONTINUOUS_LEGEND_HEIGHT,
    DEFAULT_CONTINUOUS_LEGEND_WIDTH,
} from "@/react/legend/shared/legendConstants";
export type { LegendContinuousSvgProps } from "@/react/legend/shared/legendTypes";

let nextLegendGradientId = 0;

function createLegendGradientId() {
    return `legend-gradient-${globalThis.crypto?.randomUUID?.() ?? ++nextLegendGradientId}`;
}

/**
 * Reusable continuous gradient legend SVG (bar + optional axis).
 * This is intentionally independent of color-legend wrapper logic so other
 * legend types can reuse the same renderer in future phases.
 */
export default function LegendContinuousSvg({
    label,
    colors,
    range,
    width = DEFAULT_CONTINUOUS_LEGEND_WIDTH,
    height = DEFAULT_CONTINUOUS_LEGEND_HEIGHT,
    activeRange = null,
    onRangeChange,
}: LegendContinuousSvgProps) {
    // SVG ids are document-global, and legends may be rendered through separate React roots.
    const gradientIdRef = useRef<string>(createLegendGradientId());
    const gradientId = gradientIdRef.current;
    const axisRef = useRef<SVGGElement>(null);
    const svgRef = useRef<SVGSVGElement>(null);
    const [observedWidth, setObservedWidth] = useState(width);
    const layout = getContinuousLegendLayout(observedWidth, Boolean(label));
    const colorStops = getGradientStops(colors);
    const formattedLabel = label
        ? formatLegendLabel(label, layout.labelMaxWidth)
        : null;
    const brush = useMemo(
        () =>
            onRangeChange || activeRange
                ? {
                      value: activeRange,
                      setValue: onRangeChange ?? (() => {}),
                      minMax: range,
                  }
                : undefined,
        [activeRange, onRangeChange, range],
    );

    useBrushX(svgRef, brush, "linear", {
        layout: {
            x: layout.barX,
            y: layout.barY - 3,
            width: layout.axisWidth,
            height: layout.barHeight + 6,
            handleSize: 5,
        },
        style: {
            selectionFill: "rgba(0, 0, 0, 0.18)",
            selectionStroke: "currentColor",
            handleFill: "currentColor",
        },
        updateOn: "end",
    });

    useLayoutEffect(() => {
        const container = svgRef.current?.parentElement;
        if (!container) {
            return;
        }
        const updateLayout = () => {
            setObservedWidth(container.clientWidth || width);
        };
        updateLayout();
        if (typeof ResizeObserver === "undefined") {
            return;
        }
        const observer = new ResizeObserver(updateLayout);
        observer.observe(container);
        return () => observer.disconnect();
    }, [width]);

    useLayoutEffect(() => {
        const axisG = axisRef.current;
        if (!axisG) {
            return;
        }
        const scale = scaleLinear()
            .domain([range[0], range[1]])
            .range([0, layout.axisWidth]);
        const axis = axisBottom(scale)
            .tickFormat(formatContinuousTick)
            .ticks(layout.tickCount);

        const selection = select(axisG);
        selection.call(axis);
        selection
            .selectAll("text")
            .style("text-anchor", "end")
            .attr("dx", "-.4em")
            .attr("dy", ".4em")
            .attr("transform", "rotate(-45)");

        return () => {
            selection.selectAll("*").remove();
        };
    }, [layout.axisWidth, layout.tickCount, range]);

    return (
        <svg
            ref={svgRef}
            height="100%"
            width="100%"
            className="relative overflow-hidden"
        >
            <g>
                {formattedLabel ? (
                    <g className="legend-continuous-drag-handle">
                        <rect
                            x={0}
                            y={0}
                            width="100%"
                            height={layout.barY}
                            style={{
                                fill: "transparent",
                                pointerEvents: "all",
                            }}
                        />
                        {formattedLabel.truncated ? (
                            <title>{formattedLabel.full}</title>
                        ) : null}
                        <text
                            x={10}
                            y={3}
                            aria-label={formattedLabel.full}
                            alignmentBaseline="hanging"
                            className="fill-current text-xs font-bold"
                        >
                            {formattedLabel.display}
                        </text>
                    </g>
                ) : null}
                <defs>
                    <linearGradient
                        id={gradientId}
                        x1="0%"
                        y1="0%"
                        x2="100%"
                        y2="0%"
                        spreadMethod="pad"
                    >
                        {colorStops.map((stop) => (
                            <stop
                                key={`${stop.offset}-${stop.index}`}
                                offset={stop.offset}
                                stopColor={stop.color}
                                stopOpacity={1}
                            />
                        ))}
                    </linearGradient>
                </defs>
                <rect
                    x={layout.barX}
                    y={layout.barY}
                    width={layout.axisWidth}
                    height={layout.barHeight}
                    fill={`url(#${gradientId})`}
                />
                <g
                    ref={axisRef}
                    transform={`translate(${layout.barX},${layout.barY + layout.barHeight})`}
                />
            </g>
        </svg>
    );
}
