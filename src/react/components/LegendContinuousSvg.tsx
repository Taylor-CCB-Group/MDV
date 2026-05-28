import { useId, useLayoutEffect, useRef } from "react";
import { axisBottom } from "d3-axis";
import { scaleLinear } from "d3-scale";
import { select } from "d3-selection";

export type LegendContinuousSvgProps = {
    label: string;
    colors: string[];
    range: [number, number];
    width?: number;
    height?: number;
};

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
    height = 45,
}: LegendContinuousSvgProps) {
    const gradientId = useId().replace(/:/g, "");
    const axisRef = useRef<SVGGElement>(null);
    const len = colors.length;
    const colorPct = colors.map((color, i) => {
        // Use proportional stops across [0,100] and keep keys unique via index.
        // Floor() caused repeated offsets (e.g. multiple "0%") for long palettes.
        const pct = len <= 1 ? 0 : (i / (len - 1)) * 100;
        return [`${pct}%`, color, i] as const;
    });

    let barHeight = height;
    if (label) {
        barHeight -= 10;
    }
    if (range) {
        barHeight -= 25;
    }

    useLayoutEffect(() => {
        const axisG = axisRef.current;
        if (!axisG) {
            return;
        }
        const scale = scaleLinear().domain([range[0], range[1]]).range([0, width - 20]);
        const axis = axisBottom(scale)
            .tickFormat((v) =>
                Number(v) >= 10000
                    ? Number.parseFloat(String(v)).toPrecision(2)
                    : String(v),
            )
            .ticks(Math.ceil(width / 20));

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
    }, [range, width]);

    return (
        <svg
            height={height}
            width={width}
            className="relative"
        >
            <g>
                {label ? (
                    <text
                        x={10}
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
                    y={label ? 10 : 0}
                    width={width - 20}
                    height={barHeight}
                    fill={`url(#${gradientId})`}
                />
                <g
                    ref={axisRef}
                    transform={`translate(10,${height - 25})`}
                />
            </g>
        </svg>
    );
}
