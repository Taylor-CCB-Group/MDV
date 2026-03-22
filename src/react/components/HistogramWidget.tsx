import * as d3 from "d3";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useDebounce } from "use-debounce";
import { useTheme } from "../hooks";

type Range = [number, number];
export type HistogramScaleType = "linear" | "log";

export type HistogramBrushConfig = {
    value: Range | null;
    setValue: (value: Range | null) => void;
    minMax: Range;
};

export type HistogramLayer = {
    id: string;
    data: number[];
    color: string;
    variant?: "bars" | "markers" | "line";
    inset?: number;
    widthFactor?: number;
    radius?: number;
    hidden?: boolean;
};

export type HistogramMarker = {
    id: string;
    value: number;
    color: string;
    hidden?: boolean;
};

export type HistogramScaleControls = {
    xLabel: string;
    yLabel: string;
    onToggleX: () => void;
    onToggleY: () => void;
};

type HistogramWidgetProps = {
    layers: HistogramLayer[];
    width: number;
    height: number;
    bins: number;
    binEdges?: number[];
    xScaleType?: HistogramScaleType;
    yScaleType?: HistogramScaleType;
    brush?: HistogramBrushConfig;
    markers?: HistogramMarker[];
    scaleControls?: HistogramScaleControls;
    onVisibleOnce?: () => void;
    rootMargin?: string;
};

const formatBrushValue = d3.format(".4~g");

const createScale = (
    scaleType: HistogramScaleType,
    domain: Range,
    range: Range,
) =>
    scaleType === "log"
        ? d3.scaleSymlog().domain(domain).range(range)
        : d3.scaleLinear().domain(domain).range(range);

const useBrushX = (
    ref: React.RefObject<SVGSVGElement>,
    brushConfig: HistogramBrushConfig | undefined,
    histoWidth: number,
    histoHeight: number,
    xScaleType: HistogramScaleType,
) => {
    const brushRef = useRef<ReturnType<typeof d3.brushX> | null>(null);
    const minMax = brushConfig?.minMax;
    const setValue = brushConfig?.setValue;
    const dark = useTheme() === "dark";

    useEffect(() => {
        if (!ref.current || !minMax || !setValue) return;

        const svg = d3.select(ref.current);
        const xScale = createScale(xScaleType, minMax, [0, histoWidth]);
        const brush = d3.brushX()
            .handleSize(0.5)
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
            .attr("fill", dark ? "rgba(255, 255, 255, 0.08)" : "rgba(15, 23, 42, 0.08)")
            .attr("stroke", dark ? "rgba(229, 231, 235, 0.95)" : "rgba(15, 23, 42, 0.7)")
            .attr("vector-effect", "non-scaling-stroke");
        brushGroup
            .selectAll(".handle")
            .attr("fill", dark ? "rgba(229, 231, 235, 0.92)" : "rgba(15, 23, 42, 0.88)")
            .attr("stroke", "none")
            .attr("vector-effect", "non-scaling-stroke");

        return () => {
            svg.select(".brush").remove();
        };
    }, [dark, ref, histoWidth, histoHeight, minMax, setValue, xScaleType]);

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
        const xScale = createScale(xScaleType, minMax, [0, histoWidth]);

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
    }, [histoWidth, minMax, ref, xScaleType]);

    useEffect(() => {
        setBrushValue(debouncedValue);
    }, [debouncedValue, setBrushValue]);
};

export default function HistogramWidget({
    layers,
    width,
    height,
    bins,
    binEdges,
    xScaleType = "linear",
    yScaleType = "linear",
    brush,
    markers = [],
    scaleControls,
    onVisibleOnce,
    rootMargin = "0px 0px 100px 0px",
}: HistogramWidgetProps) {
    const containerRef = useRef<HTMLDivElement>(null);
    const ref = useRef<SVGSVGElement>(null);
    const dark = useTheme() === "dark";
    const [isVisible, setIsVisible] = useState(
        typeof IntersectionObserver === "undefined",
    );
    const [hasTriggeredCallback, setHasTriggeredCallback] = useState(false);

    useEffect(() => {
        if (!containerRef.current || typeof IntersectionObserver === "undefined") return;
        const observer = new IntersectionObserver((entries) => {
            setIsVisible(entries[0]?.isIntersecting ?? false);
        });
        observer.observe(containerRef.current);
        return () => observer.disconnect();
    }, []);

    useEffect(() => {
        if (!containerRef.current || !onVisibleOnce || hasTriggeredCallback) return;
        if (typeof IntersectionObserver === "undefined") {
            setHasTriggeredCallback(true);
            onVisibleOnce();
            return;
        }
        const observer = new IntersectionObserver((entries) => {
            if (!entries[0]?.isIntersecting || hasTriggeredCallback) return;
            setHasTriggeredCallback(true);
            onVisibleOnce();
        }, { rootMargin });
        observer.observe(containerRef.current);
        return () => observer.disconnect();
    }, [onVisibleOnce, hasTriggeredCallback, rootMargin]);

    useBrushX(ref, isVisible ? brush : undefined, width, height, xScaleType);

    const padding = 2;
    const visibleLayers = useMemo(
        () => (isVisible ? layers.filter((layer) => !layer.hidden) : []),
        [isVisible, layers],
    );
    const maxValue = Math.max(1, ...visibleLayers.flatMap((layer) => layer.data));
    const xDomain = brush?.minMax ?? ([0, bins] as Range);
    const xScale = useMemo(
        () => createScale(xScaleType, xDomain, [0, width]),
        [width, xDomain, xScaleType],
    );
    const yScale = useMemo(
        () => createScale(yScaleType, [0, maxValue], [height - padding, padding]),
        [height, maxValue, padding, yScaleType],
    );

    const resolvedBinEdges = useMemo(() => {
        if (!isVisible) return [];
        if (binEdges && binEdges.length === bins + 1) {
            return binEdges;
        }
        return Array.from({ length: bins + 1 }, (_, index) =>
            xDomain[0] + ((xDomain[1] - xDomain[0]) * index) / bins,
        );
    }, [binEdges, bins, isVisible, xDomain]);

    const createBars = useCallback((data: number[]) => data.map((count, index) => {
        const startValue = resolvedBinEdges[index];
        const endValue = resolvedBinEdges[index + 1];
        const x0 = xScale(startValue);
        const x1 = xScale(endValue);
        const y = yScale(count);
        return {
            x: Math.min(x0, x1),
            y,
            width: Math.max(0.5, Math.abs(x1 - x0)),
            height: Math.max(0, height - padding - y),
        };
    }), [height, padding, resolvedBinEdges, xScale, yScale]);

    const brushValueLabel = useMemo(() => {
        if (!brush?.value) return null;
        const [start, end] = brush.value[0] <= brush.value[1]
            ? brush.value
            : [brush.value[1], brush.value[0]];
        return `${formatBrushValue(start)} - ${formatBrushValue(end)}`;
    }, [brush?.value]);

    return (
        <div className="relative w-full" ref={containerRef}>
            {brushValueLabel ? (
                <div
                    className="pointer-events-none absolute left-1 top-1 z-10 rounded px-1 py-0 text-[9px] shadow-sm"
                    style={{
                        border: `1px solid ${dark ? "rgba(71, 85, 105, 0.9)" : "rgba(203, 213, 225, 0.95)"}`,
                        background: dark ? "rgba(15, 23, 42, 0.9)" : "rgba(255, 255, 255, 0.92)",
                        color: dark ? "rgba(241, 245, 249, 0.95)" : "rgba(15, 23, 42, 0.92)",
                    }}
                >
                    {brushValueLabel}
                </div>
            ) : null}
            {scaleControls ? (
                <div className="pointer-events-none absolute right-1 top-1 z-10 flex items-center gap-1 text-[10px] opacity-75">
                    <button
                        type="button"
                        className="pointer-events-auto rounded border border-[hsl(var(--border))] bg-[hsl(var(--background)/0.9)] px-1 py-0 text-[9px] text-[hsl(var(--foreground))] shadow-sm"
                        onClick={scaleControls.onToggleX}
                    >
                        X:{scaleControls.xLabel}
                    </button>
                    <button
                        type="button"
                        className="pointer-events-auto rounded border border-[hsl(var(--border))] bg-[hsl(var(--background)/0.9)] px-1 py-0 text-[9px] text-[hsl(var(--foreground))] shadow-sm"
                        onClick={scaleControls.onToggleY}
                    >
                        Y:{scaleControls.yLabel}
                    </button>
                </div>
            ) : null}
            {isVisible ? (
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
                                        x1={bar.x + bar.width / 2}
                                        x2={bar.x + bar.width / 2}
                                        y1={height - padding - 2}
                                        y2={Math.max(padding, bar.y)}
                                        stroke={layer.color}
                                        strokeOpacity={0.2}
                                        strokeWidth={Math.max(0.7, bar.width * (layer.widthFactor ?? 0.18))}
                                        strokeDasharray="1 3"
                                        strokeLinecap="round"
                                        vectorEffect="non-scaling-stroke"
                                    />
                                );
                            });
                        }
                        if (layer.variant === "line") {
                            const points = bars
                                .map((bar, index) =>
                                    `${bar.x + bar.width / 2},${layer.data[index] === 0 ? height - padding : bar.y}`,
                                )
                                .join(" ");
                            return (
                                <polyline
                                    key={layer.id}
                                    points={points}
                                    fill="none"
                                    stroke={layer.color}
                                    strokeWidth={1.5}
                                    vectorEffect="non-scaling-stroke"
                                />
                            );
                        }
                        return bars.map((bar, index) => (
                            <rect
                                key={`${layer.id}-${index}`}
                                x={bar.x + bar.width * (layer.inset ?? 0)}
                                y={bar.y}
                                width={Math.max(0.4, bar.width * (layer.widthFactor ?? 1))}
                                height={Math.max(0, bar.height)}
                                fill={layer.color}
                                rx={layer.radius ?? 0}
                            />
                        ));
                    })}
                    {markers
                        .filter((marker) => !marker.hidden)
                        .map((marker) => {
                            const x = xScale(marker.value);
                            return (
                                <line
                                    key={marker.id}
                                    x1={x}
                                    x2={x}
                                    y1={padding + 4}
                                    y2={height - padding - 4}
                                    stroke={marker.color}
                                    strokeOpacity={0.22}
                                    strokeWidth={1}
                                    strokeLinecap="round"
                                    vectorEffect="non-scaling-stroke"
                                />
                            );
                        })}
                </svg>
            ) : (
                <div
                    aria-hidden="true"
                    style={{
                        height,
                        borderRadius: 4,
                        border: `1px solid ${dark ? "rgba(71, 85, 105, 0.45)" : "rgba(203, 213, 225, 0.9)"}`,
                        background: dark
                            ? "linear-gradient(180deg, rgba(30, 41, 59, 0.45), rgba(15, 23, 42, 0.2))"
                            : "linear-gradient(180deg, rgba(248, 250, 252, 0.95), rgba(241, 245, 249, 0.7))",
                    }}
                />
            )}
        </div>
    );
}
