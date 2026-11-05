import { type ComponentType, createElement } from "react";
import type { Root } from "react-dom/client";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import { createMdvPortal } from "@/react/react_utils";
import { makeDraggable, makeResizable } from "@/utilities/Elements.js";
import { createEl } from "@/utilities/ElementsTyped";

export type LegendWrapperPosition = {
    left: string | number;
    top: string | number;
    /** Size the user last resized the legend to, as [width, height]. */
    size?: [number, number];
};

export type LegendWrapperComponentProps<TSpec> = {
    spec: TSpec;
    /** Called after legend DOM is committed (wrapper attaches drag/resize here). */
    onLayoutReady?: () => void;
};

export type LegendWrapperInteractions = {
    /** Optional CSS selector for drag handle; defaults to full wrapper. */
    dragHandle?: string;
    /** Whether resize handles should be attached. */
    resizable?: boolean;
};

/**
 * Generic React legend wrapper for chart overlays.
 *
 * This wrapper intentionally has no color-legend-specific behavior. It can render
 * any legend component + spec pair (color, fraction, node/link, etc.) while
 * preserving shared wrapper concerns:
 * - positioned overlay container inside chart content
 * - portal lifecycle + cleanup
 * - drag/resize wiring after React commit
 * - offsetLeft/offsetTop persistence for config + SVG export
 */
export default class LegendWrapper<TSpec, T extends BaseConfig = BaseConfig> {
    private chart: BaseChart<T>;
    private wrapper: HTMLDivElement | null = null;
    private root: Root | null = null;
    private currentSpec: TSpec | null = null;
    private interactionsAttached = false;
    private interactions: LegendWrapperInteractions = {};
    private userSize: [number, number] | null = null;
    private contentObserver: ResizeObserver | null = null;

    constructor(chart: BaseChart<T>) {
        this.chart = chart;
    }

    getWrapperElement(): HTMLDivElement | null {
        return this.wrapper;
    }

    /** The size the user resized the legend to, or null if they haven't. */
    getUserSize(): [number, number] | null {
        return this.userSize;
    }

    render(
        spec: TSpec,
        position: LegendWrapperPosition,
        component: ComponentType<LegendWrapperComponentProps<TSpec>>,
        interactions: LegendWrapperInteractions = {},
    ): void {
        this.unmount();

        this.wrapper = createEl(
            "div",
            {
                classes: ["legend-container"],
                styles: {
                    position: "absolute",
                    padding: "0.4em",
                    border: "1px solid rgba(128, 128, 128, 0.65)",
                    boxShadow: "0 1px 4px rgba(0, 0, 0, 0.35), inset 0 0 0 1px rgba(255, 255, 255, 0.16)",
                },
            },
            this.chart.contentDiv,
        );

        this.applyPosition(position);
        this.currentSpec = spec;
        this.interactions = interactions;
        this.interactionsAttached = false;
        this.userSize = position.size ?? null;

        const wrapper = this.wrapper;
        this.root = createMdvPortal(
            createElement(component, {
                spec,
                onLayoutReady: () => this.attachInteractions(),
            }),
            wrapper,
            this.chart,
        );
    }

    private attachInteractions(): void {
        const wrapper = this.wrapper;
        const spec = this.currentSpec;
        if (!wrapper || !spec) {
            return;
        }
        // The legend component sets its default size in a layout effect, which can
        // run more than once (StrictMode), so the saved size and the fit-inside step
        // are reapplied every time; only the drag/resize wiring happens once.
        if (this.userSize && this.interactions.resizable) {
            wrapper.style.width = `${this.userSize[0]}px`;
            wrapper.style.height = `${this.userSize[1]}px`;
        }
        this.keepInsideChart(wrapper);
        if (this.interactionsAttached) {
            return;
        }
        const doc = this.chart.__doc__;
        makeDraggable(wrapper, {
            handle: this.interactions.dragHandle,
            contain: this.chart.contentDiv,
            doc,
        });
        if (this.interactions.resizable) {
            makeResizable(wrapper, {
                doc,
                onresizeend: (width: number, height: number) => {
                    this.userSize = [width, height];
                },
            });
        }
        // Refit when the chart changes size, including when it had no size yet above.
        if (typeof ResizeObserver !== "undefined") {
            this.contentObserver = new ResizeObserver(() => {
                if (this.wrapper) this.keepInsideChart(this.wrapper);
            });
            this.contentObserver.observe(this.chart.contentDiv);
        }
        this.interactionsAttached = true;
    }

    /** Moves the legend back inside the chart if it would hang over an edge. */
    private keepInsideChart(wrapper: HTMLDivElement): void {
        const content = this.chart.contentDiv;
        if (!content.clientWidth || !content.clientHeight) {
            return;
        }
        const maxLeft = Math.max(0, content.clientWidth - wrapper.offsetWidth);
        const maxTop = Math.max(0, content.clientHeight - wrapper.offsetHeight);
        const left = Math.min(Math.max(0, wrapper.offsetLeft), maxLeft);
        const top = Math.min(Math.max(0, wrapper.offsetTop), maxTop);
        if (left !== wrapper.offsetLeft) wrapper.style.left = `${left}px`;
        if (top !== wrapper.offsetTop) wrapper.style.top = `${top}px`;
    }

    unmount(): void {
        this.currentSpec = null;
        this.interactionsAttached = false;
        this.interactions = {};
        this.contentObserver?.disconnect();
        this.contentObserver = null;
        if (this.root) {
            this.root.unmount();
            this.root = null;
        }
        if (this.wrapper) {
            this.wrapper.remove();
            this.wrapper = null;
        }
    }

    private applyPosition(position: LegendWrapperPosition): void {
        if (!this.wrapper) {
            return;
        }
        this.wrapper.style.left = typeof position.left === "number" ? `${position.left}px` : String(position.left);
        this.wrapper.style.top = typeof position.top === "number" ? `${position.top}px` : String(position.top);
    }
}
