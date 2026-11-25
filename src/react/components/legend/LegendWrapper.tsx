import { type ComponentType, useSyncExternalStore } from "react";
import type { Root } from "react-dom/client";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import LegendFrame, {
    type LegendFrameProps,
    type LegendFrameStart,
    type LegendSummary,
} from "@/react/components/legend/LegendFrame";
import { createMdvPortal } from "@/react/react_utils";
import { createEl } from "@/utilities/ElementsTyped";

export type { LegendSummary } from "@/react/components/legend/LegendFrame";

export type LegendWrapperPosition = {
    left: string | number;
    /** Where the legend starts from the top; or give `bottom` to start from the bottom instead. */
    top?: string | number;
    bottom?: string | number;
    /** Size the user last resized the legend to, as [width, height]. */
    size?: [number, number];
};

/**
 * What a legend component is given. Where it's drawn and how much room it has come from
 * `useLegendFrame()`, and it says how big it wants to be with `useLegendPreferredSize()`.
 */
export type LegendWrapperComponentProps<TSpec> = {
    spec: TSpec;
};

export type LegendWrapperInteractions = {
    /** Optional CSS selector for drag handle; defaults to full wrapper. */
    dragHandle?: string;
    /** Whether resize handles should be attached. */
    resizable?: boolean;
    /** Shown instead of the legend while it's minimised. */
    summary?: LegendSummary;
    /** Whether this legend is minimised. Each legend is minimised on its own. */
    minimized?: boolean;
    /** Called by the legend's minimise button and by its minimised form; without it there's no button. */
    onMinimizedChange?: (minimized: boolean) => void;
    /** Called when the user has dragged the legend (or its minimised form) somewhere new. */
    onMoved?: (position: [number, number]) => void;
};

/** Where a chart wants its legends: over the plot (default) or in the strip above it. */
export type LegendPlacement = {
    onTop: boolean;
    /** Holds the legends when they're on top. */
    band: HTMLElement;
    setLayout: (change: { position?: "overlay" | "top" }) => void;
    /** Called when a legend is added to or removed from the strip. */
    onBandChanged: () => void;
};

type LegendSnapshot<TSpec> = { spec: TSpec; summary: LegendSummary };

type LegendSource<TSpec> = {
    subscribe: (listener: () => void) => () => void;
    getSnapshot: () => LegendSnapshot<TSpec> | null;
};

type FrameSettings = Omit<LegendFrameProps, "summary" | "children">;

/** The legend in its frame, showing the wrapper's current spec, so `update()` changes it without a remount. */
function LegendWrapperView<TSpec>({
    source,
    component: Legend,
    frame,
}: {
    source: LegendSource<TSpec>;
    component: ComponentType<LegendWrapperComponentProps<TSpec>>;
    frame: FrameSettings;
}) {
    const snapshot = useSyncExternalStore(source.subscribe, source.getSnapshot);
    if (snapshot === null) return null;
    return (
        <LegendFrame {...frame} summary={snapshot.summary}>
            <Legend spec={snapshot.spec} />
        </LegendFrame>
    );
}

function toPixels(value: string | number | undefined): number | undefined {
    if (value === undefined) return undefined;
    if (typeof value === "number") return value;
    const pixels = Number.parseFloat(value);
    return Number.isNaN(pixels) ? undefined : pixels;
}

/**
 * Shows a chart legend (color, fraction, contour, ...) in a `LegendFrame`, mounted from code
 * outside React: the chart says what to show and where, and the frame does the rest (dragging,
 * resizing, minimising, the strip above the plot). The legend is mounted in a host element the
 * chart can hold on to straight away, while React draws into it on its own schedule.
 */
export default class LegendWrapper<TSpec, T extends BaseConfig = BaseConfig> {
    private chart: BaseChart<T>;
    private host: HTMLDivElement | null = null;
    private root: Root | null = null;
    /** How the chart placed this legend when it was last drawn; null while it isn't shown. */
    private placement: LegendPlacement | null = null;
    private snapshot: LegendSnapshot<TSpec> | null = null;
    private listeners = new Set<() => void>();
    /** Where the legend was asked to float, until React has drawn it; null in the strip or if placed from the bottom. */
    private position: [number, number] | null = null;
    private userSize: [number, number] | null = null;
    /**
     * Whether the legend was just minimised or opened with its own button while it had keyboard
     * focus; the next time it's drawn that way, focus goes to its new button. Null otherwise.
     */
    private focusWhenMinimized: boolean | null = null;

    constructor(chart: BaseChart<T>) {
        this.chart = chart;
    }

    private subscribe = (listener: () => void) => {
        this.listeners.add(listener);
        return () => {
            this.listeners.delete(listener);
        };
    };

    private getSnapshot = () => this.snapshot;

    /** The element holding the legend, for code that needs to find, remove or export it. */
    getWrapperElement(): HTMLDivElement | null {
        return this.host;
    }

    /** The size the user resized the legend to, or null if they haven't. */
    getUserSize(): [number, number] | null {
        return this.userSize;
    }

    /**
     * Where the legend sits over the plot, to save with the view. Null when it isn't shown or
     * sits in the strip above the plot, where its offsets don't say where it floats.
     */
    getFloatingPosition(): [number, number] | null {
        if (!this.host || this.placement?.onTop) return null;
        // where the frame is now, once React has drawn it (the frame keeps it inside the chart)
        const frame = this.host.querySelector<HTMLElement>(".legend-container");
        if (frame) return [frame.offsetLeft, frame.offsetTop];
        return this.position;
    }

    render(
        spec: TSpec,
        position: LegendWrapperPosition,
        component: ComponentType<LegendWrapperComponentProps<TSpec>>,
        interactions: LegendWrapperInteractions = {},
    ): void {
        const minimized = interactions.minimized === true;
        const focusToggle = this.focusWhenMinimized === minimized;
        this.focusWhenMinimized = null;
        this.unmount();

        // charts without the method (and test stand-ins) keep legends floating over the plot
        const placement = this.chart.getLegendPlacement?.() ?? null;
        this.placement = placement;
        const inBand = placement?.onTop ?? false;

        const start: LegendFrameStart = {
            left: toPixels(position.left) ?? 0,
            top: toPixels(position.top),
            bottom: toPixels(position.bottom),
        };
        this.position = !inBand && start.top !== undefined ? [start.left, start.top] : null;
        this.userSize = position.size ?? null;
        this.snapshot = { spec, summary: interactions.summary ?? { label: "Legend" } };

        // The host has no box of its own (display: contents), so the frame inside it is laid out
        // as a child of the chart's content area, or as an item in the strip.
        const host = createEl(
            "div",
            { classes: ["legend-host"] },
            inBand && placement ? placement.band : this.chart.contentDiv,
        );
        this.host = host;
        const frame: FrameSettings = {
            inBand,
            container: this.chart.contentDiv,
            start,
            savedSize: position.size ?? null,
            resizable: interactions.resizable,
            dragHandle: interactions.dragHandle,
            minimized,
            onMinimizedChange: this.keepingFocus(interactions.onMinimizedChange),
            onPutOnTop: placement && !inBand ? () => placement.setLayout({ position: "top" }) : undefined,
            onMoved: interactions.onMoved,
            onResized: (size) => {
                // a resize finishing on a legend that's since been redrawn or removed doesn't count
                if (this.host === host) this.userSize = size;
            },
            focusToggle,
        };
        this.root = createMdvPortal(
            <LegendWrapperView
                source={{ subscribe: this.subscribe, getSnapshot: this.getSnapshot }}
                component={component}
                frame={frame}
            />,
            host,
            this.chart,
        );
        placement?.onBandChanged();
    }

    /** Notes whether the legend had keyboard focus when it's minimised or opened, before the chart redraws it. */
    private keepingFocus(
        onMinimizedChange: ((minimized: boolean) => void) | undefined,
    ): ((minimized: boolean) => void) | undefined {
        if (!onMinimizedChange) return undefined;
        return (minimized) => {
            const focused = this.host?.ownerDocument.activeElement;
            this.focusWhenMinimized = focused && this.host?.contains(focused) ? minimized : null;
            onMinimizedChange(minimized);
        };
    }

    /**
     * Shows new contents in the legend that's already drawn, without redrawing it, so its own
     * state (e.g. a hovered item) is kept. A minimised legend's button shows the new summary.
     */
    update(spec: TSpec, summary?: LegendSummary): void {
        if (!this.snapshot) return;
        this.snapshot = { spec, summary: summary ?? this.snapshot.summary };
        for (const listener of this.listeners) listener();
    }

    unmount(): void {
        this.snapshot = null;
        this.listeners.clear();
        const root = this.root;
        this.root = null;
        if (root) {
            // Legends are often redrawn from another React tree's effects (the chart's own
            // component), where unmounting a root right away makes React warn about a race.
            queueMicrotask(() => root.unmount());
        }
        this.host?.remove();
        this.host = null;
        const placement = this.placement;
        this.placement = null;
        placement?.onBandChanged();
    }
}
