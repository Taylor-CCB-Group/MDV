import * as d3 from "d3";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useTheme } from "../hooks";
import {
    createBrushScale,
    useBrushX,
    type BrushXConfig,
    type BrushXScaleType,
    type Range,
} from "@/react/components/histogram/useBrushX";

export type HistogramScaleType = BrushXScaleType;
export type HistogramBrushConfig = BrushXConfig;

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

    useBrushX(ref, isVisible ? brush : undefined, xScaleType, {
        layout: {
            width,
            height: height + 4,
            y: -2,
        },
        style: {
            selectionFill: dark
                ? "rgba(255, 255, 255, 0.08)"
                : "rgba(15, 23, 42, 0.08)",
            selectionStroke: dark
                ? "rgba(229, 231, 235, 0.95)"
                : "rgba(15, 23, 42, 0.7)",
            handleFill: dark
                ? "rgba(229, 231, 235, 0.92)"
                : "rgba(15, 23, 42, 0.88)",
        },
    });

    const padding = 2;
    const visibleLayers = useMemo(
        () => (isVisible ? layers.filter((layer) => !layer.hidden) : []),
        [isVisible, layers],
    );
    const maxValue = Math.max(1, ...visibleLayers.flatMap((layer) => layer.data));
    const xDomain = brush?.minMax ?? ([0, bins] as Range);
    const xScale = useMemo(
        () => createBrushScale(xScaleType, xDomain, [0, width]),
        [width, xDomain, xScaleType],
    );
    const yScale = useMemo(
        () =>
            createBrushScale(yScaleType, [0, maxValue], [
                height - padding,
                padding,
            ]),
        [height, maxValue, yScaleType],
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
    }), [height, resolvedBinEdges, xScale, yScale]);

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
