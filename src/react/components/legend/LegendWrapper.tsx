import { createElement, type ComponentType } from "react";
import type { Root } from "react-dom/client";
import { createMdvPortal } from "@/react/react_utils";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import { createEl } from "@/utilities/ElementsTyped";
import { makeDraggable, makeResizable } from "@/utilities/Elements.js";

export type LegendWrapperPosition = {
    left: string | number;
    top: string | number;
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
export default class LegendWrapper<
    TSpec,
    T extends BaseConfig = BaseConfig,
> {
    private chart: BaseChart<T>;
    private wrapper: HTMLDivElement | null = null;
    private root: Root | null = null;
    private currentSpec: TSpec | null = null;
    private interactionsAttached = false;
    private interactions: LegendWrapperInteractions = {};

    constructor(chart: BaseChart<T>) {
        this.chart = chart;
    }

    getWrapperElement(): HTMLDivElement | null {
        return this.wrapper;
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
                    boxShadow:
                        "0 1px 4px rgba(0, 0, 0, 0.35), inset 0 0 0 1px rgba(255, 255, 255, 0.16)",
                },
            },
            this.chart.contentDiv,
        );

        this.applyPosition(position);
        this.currentSpec = spec;
        this.interactions = interactions;
        this.interactionsAttached = false;

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
        if (!wrapper || !spec || this.interactionsAttached) {
            return;
        }
        const doc = this.chart.__doc__;
        makeDraggable(wrapper, {
            handle: this.interactions.dragHandle,
            contain: this.chart.contentDiv,
            doc,
        });
        if (this.interactions.resizable) {
            makeResizable(wrapper, { doc });
        }
        this.interactionsAttached = true;
    }

    unmount(): void {
        this.currentSpec = null;
        this.interactionsAttached = false;
        this.interactions = {};
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
        this.wrapper.style.left =
            typeof position.left === "number"
                ? `${position.left}px`
                : String(position.left);
        this.wrapper.style.top =
            typeof position.top === "number"
                ? `${position.top}px`
                : String(position.top);
    }
}
