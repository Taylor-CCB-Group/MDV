import * as d3 from "d3";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useDebounce } from "use-debounce";

type Range = [number, number];

export type HistogramBrushConfig = {
    value: Range | null;
    setValue: (value: Range | null) => void;
    minMax: Range;
};

export type HistogramLayer = {
    id: string;
    data: number[];
    color: string;
    variant?: "bars" | "markers";
    inset?: number;
    widthFactor?: number;
    radius?: number;
    hidden?: boolean;
};

type HistogramWidgetProps = {
    layers: HistogramLayer[];
    width: number;
    height: number;
    bins: number;
    brush?: HistogramBrushConfig;
    onVisibleOnce?: () => void;
    rootMargin?: string;
};

const useBrushX = (
    ref: React.RefObject<SVGSVGElement>,
    brushConfig: HistogramBrushConfig | undefined,
    histoWidth: number,
    histoHeight: number,
) => {
    const brushRef = useRef<ReturnType<typeof d3.brushX> | null>(null);
    const minMax = brushConfig?.minMax;
    const setValue = brushConfig?.setValue;

    useEffect(() => {
        if (!ref.current || !minMax || !setValue) return;

        const svg = d3.select(ref.current);
        const xScale = d3.scaleLinear().domain(minMax).range([0, histoWidth]);
        const brush = d3.brushX()
            .handleSize(1)
            .extent([[0, -2], [histoWidth, histoHeight + 2]])
            .on("brush end", (event) => {
                if (!event.sourceEvent) return;
                if (event.selection) {
                    const [start, end] = event.selection.map((x: number) =>
                        xScale.invert(x),
                    );
                    setValue([start, end]);
                } else {
                    setValue(null);
                }
            });

        brushRef.current = brush;
        const brushGroup = svg.append("g").attr("class", "brush").call(brush);
        brushGroup
            .selectAll(".selection")
            .attr("fill", "rgba(37, 99, 235, 0.12)")
            .attr("stroke", "rgba(37, 99, 235, 0.75)")
            .attr("vector-effect", "non-scaling-stroke");
        brushGroup
            .selectAll(".handle")
            .attr("fill", "rgba(37, 99, 235, 0.95)")
            .attr("vector-effect", "non-scaling-stroke");

        return () => {
            svg.select(".brush").remove();
        };
    }, [ref, histoWidth, histoHeight, minMax, setValue]);

    const [debouncedValue] = useDebounce(brushConfig?.value, 100, {
        equalityFn: (a, b) => {
            if (!a && !b) return true;
            if (!a || !b) return false;
            return a[0] === b[0] && a[1] === b[1];
        },
    });

    const setBrushValue = useCallback((value: Range | null | undefined) => {
        if (!brushRef.current || !ref.current || !minMax) return;
        const svg = d3.select(ref.current);
        const xScale = d3.scaleLinear().domain(minMax).range([0, histoWidth]);

        if (!value) {
            // @ts-ignore d3 brush typings are not worth fighting here
            svg.select(".brush").call(brushRef.current.move, null);
            return;
        }
        const [start, end] = value;
        const x0 = xScale(start);
        const x1 = xScale(end);
        // @ts-ignore d3 brush typings are not worth fighting here
        svg.select(".brush").call(brushRef.current.move, [x0, x1]);
    }, [histoWidth, minMax, ref]);

    useEffect(() => {
        setBrushValue(debouncedValue);
    }, [debouncedValue, setBrushValue]);
};

export default function HistogramWidget({
    layers,
    width,
    height,
    bins,
    brush,
    onVisibleOnce,
    rootMargin = "0px 0px 100px 0px",
}: HistogramWidgetProps) {
    const ref = useRef<SVGSVGElement>(null);
    useBrushX(ref, brush, width, height);

    const padding = 2;
    const visibleLayers = useMemo(
        () => layers.filter((layer) => !layer.hidden),
        [layers],
    );
    const maxValue = Math.max(1, ...visibleLayers.flatMap((layer) => layer.data));
    const yScale = (height - 2 * padding) / maxValue;
    const barWidth = width / bins;

    const [hasTriggeredVisible, setHasTriggeredVisible] = useState(false);
    useEffect(() => {
        if (!ref.current || !onVisibleOnce || hasTriggeredVisible) return;
        const observer = new IntersectionObserver((entries) => {
            if (!entries[0].isIntersecting || hasTriggeredVisible) return;
            setHasTriggeredVisible(true);
            onVisibleOnce();
        }, { rootMargin });
        observer.observe(ref.current);
        return () => observer.disconnect();
    }, [onVisibleOnce, hasTriggeredVisible, rootMargin]);

    const createBars = useCallback((data: number[]) => data.map((count, index) => {
        const barHeight = count * yScale;
        return {
            x: index * barWidth,
            y: height - padding - barHeight,
            height: barHeight,
        };
    }), [barWidth, height, yScale]);

    return (
        <svg
            width="100%"
            height={height}
            viewBox={`0 0 ${width} ${height}`}
            preserveAspectRatio="none"
            ref={ref}
            cursor={brush ? "move" : "default"}
        >
            {visibleLayers.map((layer) => {
                const bars = createBars(layer.data);
                if (layer.variant === "markers") {
                    return bars.map((bar, index) => {
                        if (layer.data[index] === 0) return null;
                        return (
                            <line
                                key={`${layer.id}-${index}`}
                                x1={bar.x + barWidth / 2}
                                x2={bar.x + barWidth / 2}
                                y1={height - padding}
                                y2={Math.max(padding, bar.y)}
                                stroke={layer.color}
                                strokeWidth={Math.max(1.2, barWidth * (layer.widthFactor ?? 0.35))}
                                strokeLinecap="round"
                                vectorEffect="non-scaling-stroke"
                            />
                        );
                    });
                }
                return bars.map((bar, index) => (
                    <rect
                        key={`${layer.id}-${index}`}
                        x={bar.x + barWidth * (layer.inset ?? 0)}
                        y={bar.y}
                        width={Math.max(0.4, barWidth * (layer.widthFactor ?? 1))}
                        height={Math.max(0, bar.height)}
                        fill={layer.color}
                        rx={layer.radius ?? 0}
                    />
                ));
            })}
        </svg>
    );
}
