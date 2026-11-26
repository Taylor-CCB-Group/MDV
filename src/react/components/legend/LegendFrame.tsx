import {
    type CSSProperties,
    createContext,
    type ReactNode,
    type PointerEvent as ReactPointerEvent,
    type Ref,
    type RefObject,
    useCallback,
    useContext,
    useEffect,
    useLayoutEffect,
    useMemo,
    useRef,
    useState,
} from "react";
import { flushSync } from "react-dom";
import { LEGEND_BAND_MAX_HEIGHT } from "@/react/legend/shared/legendConstants";

/** What a minimised legend shows: its title and a hint of its colours. */
export type LegendSummary = {
    label: string;
    /** A few category colours, shown as small squares. */
    swatches?: string[];
    /** How many categories the swatches leave out. */
    more?: number;
    /** Colours of a continuous scale, shown as a small bar. */
    gradient?: string[];
};

/** How big a legend's contents want to be. Without a width, the frame's width is left to the layout. */
export type LegendFrameSize = { width?: number; height: number };

/** Where a floating legend starts: from the chart's top left, or from its bottom left. */
export type LegendFrameStart = { left: number; top?: number; bottom?: number };

type LegendFrameContextValue = {
    /** True in the strip above the plot, where the legend gets an equal share of the width. */
    inBand: boolean;
    /** The width the strip gives the legend; 0 when it floats or hasn't been measured yet. */
    width: number;
    setPreferredSize: (size: LegendFrameSize | null) => void;
};

const LegendFrameContext = createContext<LegendFrameContextValue>({
    inBand: false,
    width: 0,
    setPreferredSize: () => {},
});

/** Where the legend is drawn (floating or in the strip) and the width it has there. */
export function useLegendFrame() {
    return useContext(LegendFrameContext);
}

/**
 * Tells the frame how big the legend wants to be. The frame applies it, unless the user has
 * resized the legend; in the strip it's capped in height and the width comes from the strip.
 */
export function useLegendPreferredSize(size: LegendFrameSize | null) {
    const { setPreferredSize } = useLegendFrame();
    const width = size?.width;
    const height = size?.height;
    useLayoutEffect(() => {
        setPreferredSize(height === undefined ? null : { width, height });
    }, [setPreferredSize, width, height]);
}

export type LegendFrameProps = {
    inBand: boolean;
    /** The chart's content area: a floating legend is positioned in it and kept inside it. */
    container: HTMLElement;
    start: LegendFrameStart;
    /** The size the user last resized the legend to, used only if it's resizable. */
    savedSize?: [number, number] | null;
    resizable?: boolean;
    /** CSS selector for the part the legend is dragged by; the whole legend if not given. */
    dragHandle?: string;
    minimized: boolean;
    summary: LegendSummary;
    /** Without it, the legend can't be minimised. */
    onMinimizedChange?: (minimized: boolean) => void;
    /** Moves the chart's legends into the strip above the plot; without it there's no button. */
    onPutOnTop?: () => void;
    /** Called when the user has dragged the legend somewhere new. */
    onMoved?: (position: [number, number]) => void;
    onResized?: (size: [number, number]) => void;
    /**
     * Moves keyboard focus to the minimise or + button once drawn, for a legend that was just
     * minimised or opened with its own button, so focus isn't lost when it's redrawn.
     */
    focusToggle?: boolean;
    children: ReactNode;
};

type Point = { left: number; top: number };

type ResizeDirection = "n" | "s" | "e" | "w" | "ne" | "nw" | "se" | "sw";
const RESIZE_DIRECTIONS: ResizeDirection[] = ["n", "s", "e", "w", "ne", "nw", "se", "sw"];
const MIN_RESIZED_SIZE = 40;

function clamp(value: number, min: number, max: number) {
    return Math.min(Math.max(value, min), Math.max(min, max));
}

/**
 * Follows the pointer that pressed the legend, from the document it's in (charts can be popped
 * out into another window), until it's released; other pointers are ignored. Returns a function
 * that stops following it without calling `onEnd`.
 */
function trackPointer(
    event: ReactPointerEvent<HTMLElement>,
    onMove: (dx: number, dy: number) => void,
    onEnd: () => void,
): () => void {
    const doc = event.currentTarget.ownerDocument;
    const { pointerId } = event;
    const startX = event.clientX;
    const startY = event.clientY;
    const move = (e: PointerEvent) => {
        if (e.pointerId !== pointerId) return;
        e.preventDefault();
        onMove(e.clientX - startX, e.clientY - startY);
    };
    const stop = () => {
        doc.removeEventListener("pointermove", move);
        doc.removeEventListener("pointerup", end);
        doc.removeEventListener("pointercancel", end);
    };
    const end = (e: PointerEvent) => {
        if (e.pointerId !== pointerId) return;
        stop();
        onEnd();
    };
    doc.addEventListener("pointermove", move);
    doc.addEventListener("pointerup", end);
    doc.addEventListener("pointercancel", end);
    return stop;
}

/**
 * Whether the event target is an element. Not `instanceof Element`: a chart popped out into
 * another window holds elements from both windows, and each window has its own `Element`.
 */
function isElement(target: EventTarget | null): target is Element {
    return target !== null && "closest" in target;
}

/** The element's width, followed as it changes; 0 while not `enabled`. */
function useElementWidth(ref: RefObject<HTMLElement | null>, enabled: boolean) {
    const [width, setWidth] = useState(0);
    useLayoutEffect(() => {
        const el = ref.current;
        if (!enabled || !el) {
            setWidth(0);
            return;
        }
        const measure = () => setWidth(el.clientWidth);
        measure();
        if (typeof ResizeObserver === "undefined") return;
        const observer = new ResizeObserver(measure);
        observer.observe(el);
        return () => observer.disconnect();
    }, [ref, enabled]);
    return width;
}

/**
 * A floating legend's position in the chart. It starts where it was saved, and is kept inside the
 * chart: after every render, since the legend's size can change, and whenever the chart is resized.
 * One placed from the bottom stays pinned there (so it doesn't matter how tall it turns out to be)
 * until it's dragged, resized, or has to be moved to fit.
 */
function useFloatingPosition(
    ref: RefObject<HTMLElement | null>,
    container: HTMLElement,
    start: LegendFrameStart,
    enabled: boolean,
) {
    const [position, setPosition] = useState<Point | null>(
        start.top === undefined ? null : { left: start.left, top: start.top },
    );

    const keepInside = useCallback(() => {
        const el = ref.current;
        const width = container.clientWidth;
        const height = container.clientHeight;
        if (!el || !width || !height) return;
        setPosition((current) => {
            const base = current ?? { left: el.offsetLeft, top: el.offsetTop };
            const left = clamp(base.left, 0, width - el.offsetWidth);
            const top = clamp(base.top, 0, height - el.offsetHeight);
            if (left === base.left && top === base.top) return current;
            return { left, top };
        });
    }, [ref, container]);

    useLayoutEffect(() => {
        if (enabled) keepInside();
    });

    useEffect(() => {
        if (!enabled || typeof ResizeObserver === "undefined") return;
        const observer = new ResizeObserver(keepInside);
        observer.observe(container);
        return () => observer.disconnect();
    }, [enabled, container, keepInside]);

    /** Moves the legend, keeping it inside the chart, and says where it ended up. */
    const moveTo = useCallback(
        (left: number, top: number): Point => {
            const el = ref.current;
            const moved = el
                ? {
                      left: clamp(left, 0, container.clientWidth - el.offsetWidth),
                      top: clamp(top, 0, container.clientHeight - el.offsetHeight),
                  }
                : { left, top };
            setPosition(moved);
            return moved;
        },
        [ref, container],
    );

    const style: CSSProperties = position
        ? { position: "absolute", left: position.left, top: position.top }
        : { position: "absolute", left: start.left, bottom: start.bottom ?? 0 };
    return { style, moveTo, setPosition };
}

/** A small icon button on a legend. */
export function LegendControlButton({
    icon,
    label,
    onClick,
    className,
    ref,
}: {
    icon: string;
    label: string;
    onClick: () => void;
    className?: string;
    ref?: Ref<HTMLButtonElement>;
}) {
    return (
        <button
            ref={ref}
            type="button"
            className={className ? `legend-control-btn ${className}` : "legend-control-btn"}
            title={label}
            aria-label={label}
            onClick={onClick}
        >
            <i className={icon} aria-hidden="true" />
        </button>
    );
}

/** The minimised legend's title and colours, with a + button that opens the legend again. */
function LegendSummaryPill({
    summary,
    onOpen,
    openRef,
}: {
    summary: LegendSummary;
    onOpen: () => void;
    openRef: Ref<HTMLButtonElement>;
}) {
    return (
        <div className="legend-summary" title={summary.label}>
            <span className="legend-summary-label">{summary.label}</span>
            {/* the colours give way to the + on a narrow pill, and are cut off */}
            <span className="legend-summary-previews">
                {(summary.swatches ?? []).map((color, i) => (
                    // biome-ignore lint/suspicious/noArrayIndexKey: swatches are positional and colours can repeat
                    <span key={i} className="legend-summary-swatch" style={{ backgroundColor: color }} />
                ))}
                {summary.more ? <span className="legend-summary-more">+{summary.more}</span> : null}
                {summary.gradient && summary.gradient.length > 1 ? (
                    <span
                        className="legend-summary-gradient"
                        style={{ backgroundImage: `linear-gradient(90deg, ${summary.gradient.join(", ")})` }}
                    />
                ) : null}
            </span>
            <LegendControlButton
                icon="fas fa-plus"
                label={`Show legend: ${summary.label}`}
                className="legend-summary-open"
                onClick={onOpen}
                ref={openRef}
            />
        </div>
    );
}

/**
 * The frame every chart legend is drawn in. Floating over the plot it can be dragged, resized
 * and minimised, and is kept inside the chart; in the strip above the plot it takes an equal share
 * of the width. Minimised, it's a small button with the legend's title and colours.
 */
export default function LegendFrame({
    inBand,
    container,
    start,
    savedSize = null,
    resizable = false,
    dragHandle,
    minimized: minimizedProp,
    summary,
    onMinimizedChange,
    onPutOnTop,
    onMoved,
    onResized,
    focusToggle = false,
    children,
}: LegendFrameProps) {
    const ref = useRef<HTMLDivElement>(null);
    /** The minimise button, or the + while minimised. */
    const toggleRef = useRef<HTMLButtonElement>(null);
    /** Stops the drag or resize in progress, if any. */
    const gesture = useRef<(() => void) | null>(null);
    const minimized = minimizedProp && onMinimizedChange !== undefined;
    const floating = !inBand;
    const canResize = resizable && floating && !minimized;

    const [preferred, setPreferred] = useState<LegendFrameSize | null>(null);
    const [userSize, setUserSize] = useState<[number, number] | null>(resizable ? savedSize : null);
    const bandWidth = useElementWidth(ref, inBand && !minimized);
    const { style: positionStyle, moveTo, setPosition } = useFloatingPosition(ref, container, start, floating);

    useLayoutEffect(() => {
        if (focusToggle) toggleRef.current?.focus();
    }, [focusToggle]);

    // a drag or resize still going when the legend is redrawn is dropped, not reported
    useEffect(
        () => () => {
            gesture.current?.();
            gesture.current = null;
        },
        [],
    );

    const track = (
        event: ReactPointerEvent<HTMLDivElement>,
        onMove: (dx: number, dy: number) => void,
        onEnd: () => void,
    ) => {
        // the first pointer keeps control until it's released
        if (gesture.current) return;
        gesture.current = trackPointer(event, onMove, () => {
            gesture.current = null;
            onEnd();
        });
    };

    /**
     * Where the frame is shown once `update` is drawn, along with anything that keeps it inside
     * the chart afterwards: what a drag or resize reports when it ends, so the saved position
     * matches the legend on screen.
     */
    const shownAfter = (update: () => void): [number, number] | null => {
        flushSync(update);
        const el = ref.current;
        return el ? [el.offsetLeft, el.offsetTop] : null;
    };

    const onPointerDown = (event: ReactPointerEvent<HTMLDivElement>) => {
        const el = ref.current;
        const target = event.target;
        if (!floating || !el || event.button !== 0 || !isElement(target)) return;
        if (target.closest("button, [data-resize]")) return;
        if (!minimized && dragHandle && !target.closest(dragHandle)) return;
        event.preventDefault();
        const from = { left: el.offsetLeft, top: el.offsetTop };
        let last: Point | null = null;
        track(
            event,
            (dx, dy) => {
                last = moveTo(from.left + dx, from.top + dy);
            },
            () => {
                const end = last;
                if (!end) return;
                const shown = shownAfter(() => setPosition(end));
                if (shown) onMoved?.(shown);
            },
        );
    };

    const onResizeStart = (direction: ResizeDirection) => (event: ReactPointerEvent<HTMLDivElement>) => {
        const el = ref.current;
        if (!el || event.button !== 0) return;
        event.preventDefault();
        event.stopPropagation();
        const from = { left: el.offsetLeft, top: el.offsetTop, width: el.offsetWidth, height: el.offsetHeight };
        let size: [number, number] | null = null;
        let moved: Point | null = null;
        track(
            event,
            (dx, dy) => {
                let { left, top, width, height } = from;
                if (direction.includes("e")) width = from.width + dx;
                if (direction.includes("s")) height = from.height + dy;
                if (direction.includes("w")) {
                    width = Math.max(MIN_RESIZED_SIZE, from.width - dx);
                    left = from.left + from.width - width;
                }
                if (direction.includes("n")) {
                    height = Math.max(MIN_RESIZED_SIZE, from.height - dy);
                    top = from.top + from.height - height;
                }
                size = [Math.max(MIN_RESIZED_SIZE, width), Math.max(MIN_RESIZED_SIZE, height)];
                setUserSize(size);
                if (left !== from.left || top !== from.top) moved = { left, top };
                setPosition({ left, top });
            },
            () => {
                const endSize = size;
                const endPosition = moved;
                if (endSize) onResized?.(endSize);
                if (!endPosition) return;
                const shown = shownAfter(() => {
                    if (endSize) setUserSize(endSize);
                    setPosition(endPosition);
                });
                if (shown) onMoved?.(shown);
            },
        );
    };

    const context = useMemo<LegendFrameContextValue>(
        () => ({ inBand, width: bandWidth, setPreferredSize: setPreferred }),
        [inBand, bandWidth],
    );

    if (minimized) {
        return (
            <div
                ref={ref}
                className={
                    inBand ? "legend-container legend-in-band legend-minimized" : "legend-container legend-minimized"
                }
                style={floating ? { ...positionStyle, touchAction: "none" } : undefined}
                onPointerDown={floating ? onPointerDown : undefined}
            >
                <LegendSummaryPill summary={summary} onOpen={() => onMinimizedChange?.(false)} openRef={toggleRef} />
            </div>
        );
    }

    // the size the user gave a floating legend is kept for when it floats again, not used in the strip
    const width = floating ? (userSize?.[0] ?? preferred?.width) : undefined;
    const height = floating ? (userSize?.[1] ?? preferred?.height) : preferred?.height;
    const sizeStyle: CSSProperties = {
        width,
        height: inBand && height !== undefined ? Math.min(height, LEGEND_BAND_MAX_HEIGHT) : height,
        // dragged from anywhere: touch moves it rather than scrolling the page
        touchAction: floating && !dragHandle ? "none" : undefined,
    };

    return (
        <div
            ref={ref}
            className={inBand ? "legend-container legend-in-band" : "legend-container legend-floating"}
            style={floating ? { ...positionStyle, ...sizeStyle } : sizeStyle}
            onPointerDown={floating ? onPointerDown : undefined}
        >
            <LegendFrameContext.Provider value={context}>{children}</LegendFrameContext.Provider>
            {onPutOnTop || onMinimizedChange ? (
                <div className="legend-controls">
                    {floating && onPutOnTop ? (
                        <LegendControlButton
                            icon="fas fa-window-maximize"
                            label="Put legends above the plot"
                            onClick={onPutOnTop}
                        />
                    ) : null}
                    {onMinimizedChange ? (
                        <LegendControlButton
                            icon="fas fa-minus"
                            label="Minimise legend"
                            onClick={() => onMinimizedChange(true)}
                            ref={toggleRef}
                        />
                    ) : null}
                </div>
            ) : null}
            {canResize
                ? RESIZE_DIRECTIONS.map((direction) => (
                      <div
                          key={direction}
                          className={`resizer-${direction}`}
                          data-resize={direction}
                          style={{ touchAction: "none" }}
                          onPointerDown={onResizeStart(direction)}
                      />
                  ))
                : null}
        </div>
    );
}
