import { useLayoutEffect, useRef, useState } from "react";
import { axisBottom } from "d3";
import { scaleLinear } from "d3-scale";
import { select } from "d3-selection";

let nextLegendGradientId = 0;

function createLegendGradientId() {
    return `legend-gradient-${globalThis.crypto?.randomUUID?.() ?? ++nextLegendGradientId}`;
}

export type LegendContinuousSvgProps = {
    label: string;
    colors: string[];
    range: [number, number];
    width?: number;
    height?: number;
};

export const DEFAULT_CONTINUOUS_LEGEND_HEIGHT = 65;

/**
 * Reusable continuous gradient legend SVG (bar + optional axis).
 * This is intentionally independent of color-legend wrapper logic so other
 * legend types can reuse the same renderer in future phases.
 */
export default function LegendContinuousSvg({
    label,
    colors,
    range,
    width = 120,
    height = DEFAULT_CONTINUOUS_LEGEND_HEIGHT,
}: LegendContinuousSvgProps) {
    // SVG ids are document-global, and legends may be rendered through separate React roots.
    const gradientIdRef = useRef<string>(createLegendGradientId());
    const gradientId = gradientIdRef.current;
    const axisRef = useRef<SVGGElement>(null);
    const svgRef = useRef<SVGSVGElement>(null);
    const [observedWidth, setObservedWidth] = useState(width);
    const layoutWidth = Math.max(observedWidth, 40);
    const len = colors.length;
    const colorPct = colors.map((color, i) => {
        // Use proportional stops across [0,100] and keep keys unique via index.
        // Floor() caused repeated offsets (e.g. multiple "0%") for long palettes.
        const pct = len <= 1 ? 0 : (i / (len - 1)) * 100;
        return [`${pct}%`, color, i] as const;
    });

    const barY = label ? 15 : 0;
    const barHeight = 10;

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
    }, [height, width]);

    useLayoutEffect(() => {
        const axisG = axisRef.current;
        if (!axisG) {
            return;
        }
        const scale = scaleLinear()
            .domain([range[0], range[1]])
            .range([0, layoutWidth - 20]);
        const axis = axisBottom(scale)
            .tickFormat((v) =>
                Number(v) >= 10000
                    ? Number.parseFloat(String(v)).toPrecision(2)
                    : String(v),
            )
            .ticks(Math.max(2, Math.ceil(layoutWidth / 35)));

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
    }, [layoutWidth, range]);

    return (
        <svg
            ref={svgRef}
            height="100%"
            width="100%"
            className="relative overflow-hidden"
        >
            <g>
                {label ? (
                    <text
                        x={10}
                        aria-label={label}
                        alignmentBaseline="hanging"
                        className="fill-current text-xs"
                    >
                        {label}
                    </text>
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
                        {colorPct.map(([offset, stopColor, index]) => (
                            <stop
                                key={`${offset}-${index}`}
                                offset={offset}
                                stopColor={stopColor}
                                stopOpacity={1}
                            />
                        ))}
                    </linearGradient>
                </defs>
                <rect
                    x={10}
                    y={barY}
                    width={layoutWidth - 20}
                    height={barHeight}
                    fill={`url(#${gradientId})`}
                />
                <g
                    ref={axisRef}
                    transform={`translate(10,${barY + barHeight})`}
                />
            </g>
        </svg>
    );
}
