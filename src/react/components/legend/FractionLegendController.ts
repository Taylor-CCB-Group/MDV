import { runInAction } from "mobx";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import FractionLegend from "@/react/components/legend/FractionLegend";
import LegendWrapper from "@/react/components/legend/LegendWrapper";
import type { FractionLegendSpec } from "@/react/legend/fraction_legend/types";

/** A dot plot's fraction legend settings, saved with the view. */
export type FractionLegendSettings = {
    display: boolean;
    /** Where it floats over the plot, [left, top] in pixels. */
    position?: [number, number];
    /** Older views saved the position here. */
    pos?: [number, number];
    minimized?: boolean;
};

/** Builds what the legend shows, given where it floats; null when there's nothing to show yet. */
export type FractionLegendBuilder = (position: [number, number]) => FractionLegendSpec | null;

type ShowRequest = { build: FractionLegendBuilder | null; defaultPosition: [number, number] };

/**
 * A dot plot's "% expressing" legend, shared by the original dot plot and the React one. It
 * saves where the legend floats and whether it's minimised with the view, and redraws itself when
 * the chart's legends move above or over the plot.
 */
export default class FractionLegendController<T extends BaseConfig = BaseConfig> {
    private wrapper: LegendWrapper<FractionLegendSpec, T>;
    private request: ShowRequest | null = null;
    private stopWatchingLayout: () => void;

    constructor(
        chart: BaseChart<T>,
        private settings: () => FractionLegendSettings,
    ) {
        this.wrapper = new LegendWrapper<FractionLegendSpec, T>(chart);
        this.stopWatchingLayout = chart.subscribeLegendLayout(() => this.redraw());
    }

    /**
     * Shows the legend, or hides it if it's switched off or `build` is null. It starts at
     * `defaultPosition` until the user moves it.
     */
    show(build: FractionLegendBuilder | null, defaultPosition: [number, number]): void {
        this.request = { build, defaultPosition };
        this.redraw();
    }

    /** Saves where the legend floats over the plot, e.g. before the view is saved. */
    savePosition(): void {
        const floating = this.wrapper.getFloatingPosition();
        if (floating) this.setPosition(floating);
    }

    /** Removes the legend for good, when the chart is removed. */
    remove(): void {
        this.stopWatchingLayout();
        this.request = null;
        this.wrapper.unmount();
    }

    private redraw(): void {
        this.savePosition();
        const settings = this.settings();
        const build = this.request?.build;
        if (!this.request || !build || !settings.display) {
            this.wrapper.unmount();
            return;
        }
        const position = settings.position ?? settings.pos ?? this.request.defaultPosition;
        // kept, so the legend stays put if the chart's default place for it changes
        this.setPosition(position);
        const spec = build(position);
        if (!spec) {
            this.wrapper.unmount();
            return;
        }
        // the size the user resized it to, kept across redraws
        const size = this.wrapper.getUserSize() ?? undefined;
        this.wrapper.render(spec, { left: position[0], top: position[1], size }, FractionLegend, {
            dragHandle: ".legend-drag-handle",
            resizable: true,
            summary: { label: spec.label },
            minimized: settings.minimized === true,
            onMinimizedChange: (minimized) => {
                runInAction(() => {
                    this.settings().minimized = minimized;
                });
                this.redraw();
            },
        });
    }

    private setPosition(position: [number, number]): void {
        runInAction(() => {
            this.settings().position = [position[0], position[1]];
        });
    }
}
