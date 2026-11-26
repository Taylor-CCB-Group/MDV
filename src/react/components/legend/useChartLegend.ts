import { runInAction } from "mobx";
import { type ComponentType, useEffect, useRef, useSyncExternalStore } from "react";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import LegendWrapper, { type LegendSummary, type LegendWrapperComponentProps } from "./LegendWrapper";

/** A legend's saved settings in the chart config. */
export type ChartLegendConfig = {
    display: boolean;
    /** Where the user last left it over the plot, [left, top] in pixels. */
    pos?: [number, number];
    minimized?: boolean;
};

type UseChartLegendOptions<TSpec, T extends BaseConfig> = {
    chart: BaseChart<T>;
    /** The legend's settings in the (observable) chart config. */
    legend: ChartLegendConfig;
    /** What to show; null hides the legend. Keep it stable (e.g. memoised) between renders. */
    spec: TSpec | null;
    component: ComponentType<LegendWrapperComponentProps<TSpec>>;
    summary: LegendSummary;
    /** Where it starts before the user moves it: [left, bottom] in pixels from the plot's corner. */
    defaultBottomLeft: [number, number];
};

/**
 * Shows a legend drawn by a React chart in the shared `LegendWrapper`, so it can be dragged,
 * minimised and moved above the plot like the colour legend. Its contents update in place;
 * it's only redrawn when where or how it's shown changes.
 */
export function useChartLegend<TSpec, T extends BaseConfig>({
    chart,
    legend,
    spec,
    component,
    summary,
    defaultBottomLeft,
}: UseChartLegendOptions<TSpec, T>) {
    const wrapperRef = useRef<LegendWrapper<TSpec, T> | null>(null);
    if (!wrapperRef.current) wrapperRef.current = new LegendWrapper<TSpec, T>(chart);
    // the redraw effect reads these without depending on them; kept up to date before it runs
    const latest = useRef({ spec, summary, legend, defaultBottomLeft });
    useEffect(() => {
        latest.current = { spec, summary, legend, defaultBottomLeft };
    });

    const visible = legend.display && spec !== null;
    const onTop = useSyncExternalStore(chart.subscribeLegendLayout, () => chart.legendsOnTop());
    const minimized = legend.minimized === true;

    // biome-ignore lint/correctness/useExhaustiveDependencies: redraw only when placement changes; contents update below
    useEffect(() => {
        const wrapper = wrapperRef.current;
        const { spec: current, summary: currentSummary, legend: settings, defaultBottomLeft: start } = latest.current;
        if (!wrapper || !visible || current === null) return;
        const pos = settings.pos;
        wrapper.render(current, pos ? { left: pos[0], top: pos[1] } : { left: start[0], bottom: start[1] }, component, {
            dragHandle: ".legend-drag-handle",
            summary: currentSummary,
            minimized,
            onMinimizedChange: (value) =>
                runInAction(() => {
                    settings.minimized = value;
                }),
            onMoved: (position) =>
                runInAction(() => {
                    settings.pos = position;
                }),
        });
        return () => wrapper.unmount();
    }, [visible, onTop, minimized, component]);

    useEffect(() => {
        if (spec !== null) wrapperRef.current?.update(spec, summary);
    }, [spec, summary]);
}
