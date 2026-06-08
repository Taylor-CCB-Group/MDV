import {
    useEffect,
    useLayoutEffect,
    useRef,
    useState,
    type MouseEvent,
} from "react";
import { axisBottom } from "d3";
import { scaleLinear } from "d3-scale";
import { select } from "d3-selection";
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
    const barRef = useRef<SVGRectElement>(null);
    const svgRef = useRef<SVGSVGElement>(null);
    const [observedWidth, setObservedWidth] = useState(width);
    const dragStartRef = useRef<number | null>(null);
    const cleanupDocumentRangeListenersRef = useRef<(() => void) | null>(null);
    const [draftRange, setDraftRange] = useState<[number, number] | null>(null);
    const layout = getContinuousLegendLayout(observedWidth, Boolean(label));
    const colorStops = getGradientStops(colors);
    const formattedLabel = label
        ? formatLegendLabel(label, layout.labelMaxWidth)
        : null;
    const interactive = Boolean(onRangeChange);

    const valueToX = (value: number) => {
        const span = range[1] - range[0];
        if (span === 0) {
            return layout.barX;
        }
        const fraction = Math.min(
            1,
            Math.max(0, (value - range[0]) / span),
        );
        return layout.barX + fraction * layout.axisWidth;
    };

    const getClientValue = (clientX: number) => {
        const rect = barRef.current?.getBoundingClientRect();
        if (!rect || rect.width === 0) {
            return range[0];
        }
        const fraction = Math.min(
            1,
            Math.max(0, (clientX - rect.left) / rect.width),
        );
        return range[0] + fraction * (range[1] - range[0]);
    };

    const normalizeRange = (start: number, end: number): [number, number] =>
        start <= end ? [start, end] : [end, start];

    const cleanupDocumentRangeListeners = () => {
        cleanupDocumentRangeListenersRef.current?.();
        cleanupDocumentRangeListenersRef.current = null;
    };

    const finishRangeSelection = (clientX: number) => {
        if (!onRangeChange || dragStartRef.current === null) {
            return;
        }
        const start = dragStartRef.current;
        const end = getClientValue(clientX);
        const span = Math.abs(range[1] - range[0]);
        const nextRange = normalizeRange(start, end);
        dragStartRef.current = null;
        setDraftRange(null);

        if (Math.abs(end - start) < span * 0.005) {
            if (activeRange) {
                onRangeChange(null);
            }
            return;
        }
        onRangeChange(nextRange);
    };

    const handleRangeMouseDown = (event: MouseEvent<SVGElement>) => {
        if (!onRangeChange) {
            return;
        }
        event.preventDefault();
        event.stopPropagation();
        cleanupDocumentRangeListeners();
        const value = getClientValue(event.clientX);
        dragStartRef.current = value;
        setDraftRange([value, value]);

        const ownerDocument = svgRef.current?.ownerDocument;
        if (!ownerDocument) {
            return;
        }
        const handleDocumentMouseMove = (documentEvent: globalThis.MouseEvent) => {
            if (dragStartRef.current === null) {
                return;
            }
            documentEvent.preventDefault();
            documentEvent.stopPropagation();
            setDraftRange(
                normalizeRange(
                    dragStartRef.current,
                    getClientValue(documentEvent.clientX),
                ),
            );
        };
        const handleDocumentMouseUp = (documentEvent: globalThis.MouseEvent) => {
            documentEvent.preventDefault();
            documentEvent.stopPropagation();
            cleanupDocumentRangeListeners();
            finishRangeSelection(documentEvent.clientX);
        };
        ownerDocument.addEventListener("mousemove", handleDocumentMouseMove);
        ownerDocument.addEventListener("mouseup", handleDocumentMouseUp);
        cleanupDocumentRangeListenersRef.current = () => {
            ownerDocument.removeEventListener("mousemove", handleDocumentMouseMove);
            ownerDocument.removeEventListener("mouseup", handleDocumentMouseUp);
        };
    };

    const handleRangeMouseMove = (event: MouseEvent<SVGElement>) => {
        if (dragStartRef.current === null) {
            return;
        }
        event.stopPropagation();
        setDraftRange(normalizeRange(dragStartRef.current, getClientValue(event.clientX)));
    };

    const handleRangeMouseUp = (event: MouseEvent<SVGElement>) => {
        if (!onRangeChange || dragStartRef.current === null) {
            return;
        }
        event.stopPropagation();
        cleanupDocumentRangeListeners();
        finishRangeSelection(event.clientX);
    };

    const selectedRange = draftRange ?? activeRange;
    const selectedBounds = selectedRange
        ? {
              left: valueToX(selectedRange[0]),
              right: valueToX(selectedRange[1]),
          }
        : null;

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

    useEffect(
        () => () => {
            cleanupDocumentRangeListenersRef.current?.();
            cleanupDocumentRangeListenersRef.current = null;
        },
        [],
    );

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
                        {formattedLabel.truncated ? (
                            <title>{formattedLabel.full}</title>
                        ) : null}
                        <text
                            x={10}
                            aria-label={formattedLabel.full}
                            alignmentBaseline="hanging"
                            className="fill-current text-xs"
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
                    ref={barRef}
                    x={layout.barX}
                    y={layout.barY}
                    width={layout.axisWidth}
                    height={layout.barHeight}
                    fill={`url(#${gradientId})`}
                />
                {selectedBounds ? (
                    <g pointerEvents="none">
                        <rect
                            x={layout.barX}
                            y={layout.barY}
                            width={Math.max(0, selectedBounds.left - layout.barX)}
                            height={layout.barHeight}
                            fill="rgba(0,0,0,0.45)"
                        />
                        <rect
                            x={selectedBounds.right}
                            y={layout.barY}
                            width={Math.max(
                                0,
                                layout.barX +
                                    layout.axisWidth -
                                    selectedBounds.right,
                            )}
                            height={layout.barHeight}
                            fill="rgba(0,0,0,0.45)"
                        />
                        <rect
                            x={selectedBounds.left}
                            y={layout.barY - 1}
                            width={Math.max(
                                1,
                                selectedBounds.right - selectedBounds.left,
                            )}
                            height={layout.barHeight + 2}
                            fill="none"
                            stroke="currentColor"
                            strokeWidth={1.5}
                        />
                    </g>
                ) : null}
                <g
                    ref={axisRef}
                    transform={`translate(${layout.barX},${layout.barY + layout.barHeight})`}
                />
                <rect
                    x={Math.max(0, layout.barX - 8)}
                    y={layout.barY - 4}
                    width={layout.axisWidth + 16}
                    height={layout.barHeight + 8}
                    fill="transparent"
                    role={interactive ? "slider" : undefined}
                    aria-label={
                        interactive
                            ? `Filter ${label || "continuous legend"} range`
                            : undefined
                    }
                    className={interactive ? "cursor-crosshair" : undefined}
                    onMouseDown={handleRangeMouseDown}
                    onMouseMove={handleRangeMouseMove}
                    onMouseUp={handleRangeMouseUp}
                />
            </g>
        </svg>
    );
}
