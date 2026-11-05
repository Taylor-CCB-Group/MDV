import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import BaseChart from "@/charts/BaseChart";
import type { LegendWrapperComponentProps } from "@/react/components/legend/LegendWrapper";
import type {
    ColorLegendCategoricalSpec,
    ColorLegendContinuousSpec,
    ColorLegendSpec,
} from "@/react/legend/color_legend/types";

const categorical: ColorLegendCategoricalSpec = {
    kind: "categorical",
    label: "Cell type",
    column: "cell_type",
    items: [{ color: "#ff0000", name: "T cell", value: "T cell" }],
};

const continuous: ColorLegendContinuousSpec = {
    kind: "continuous",
    label: "CD14",
    column: "cd14",
    colors: ["#000000", "#ffffff"],
    range: [0, 10],
};

type FakeOptions = {
    spec: ColorLegendSpec;
    colorLegend?: Record<string, unknown>;
    /** A legend already on the chart, as when the legend is redrawn. */
    existingLegend?: HTMLDivElement;
    userSize?: [number, number] | null;
};

/** A chart with just the parts `setColorLegend` and `getConfig` use. */
function fakeChart({ spec, colorLegend = { display: true }, existingLegend, userSize = null }: FakeOptions) {
    const legendEl = document.createElement("div");
    const wrapper = {
        render: vi.fn(),
        unmount: vi.fn(),
        getWrapperElement: vi.fn(() => legendEl),
        getUserSize: vi.fn(() => userSize),
    };
    const getDimension = vi.fn();
    const chart = Object.assign(Object.create(BaseChart.prototype), {
        config: { color_legend: colorLegend },
        legend: existingLegend,
        colorLegendWrapper: wrapper,
        colorLegendFilter: null,
        colorLegendFilterDimension: null,
        colorLegendFilterDimensionKind: null,
        dataStore: { getDimension },
        width: 400,
        height: 300,
        __doc__: document,
        getColorLegendSpec: () => spec,
        updateResetButtonVisibility: vi.fn(),
        activeQueries: { serialiseChartConfig: vi.fn(() => chart.config) },
    });
    return { chart, wrapper, getDimension };
}

/** What `setColorLegend` handed to the legend wrapper. */
function renderedLegend(wrapper: { render: ReturnType<typeof vi.fn> }) {
    const [spec, position, Component, interactions] = wrapper.render.mock.calls[0];
    const props: LegendWrapperComponentProps<ColorLegendSpec> = { spec };
    return { position, interactions, legendProps: Component(props).props };
}

beforeEach(() => {
    vi.spyOn(HTMLElement.prototype, "offsetLeft", "get").mockImplementation(function (this: HTMLElement) {
        return Number.parseFloat(this.style.left) || 0;
    });
    vi.spyOn(HTMLElement.prototype, "offsetTop", "get").mockImplementation(function (this: HTMLElement) {
        return Number.parseFloat(this.style.top) || 0;
    });
});

afterEach(() => {
    vi.restoreAllMocks();
});

describe("BaseChart color legend size", () => {
    test("passes a saved size to the legend when the chart loads", () => {
        const { chart, wrapper } = fakeChart({
            spec: categorical,
            colorLegend: { display: true, pos: [10, 20], size: [120, 90] },
        });
        chart.setColorLegend();
        expect(renderedLegend(wrapper).position).toEqual({ left: "10px", top: "20px", size: [120, 90] });
    });

    test("ignores a saved size that isn't two positive numbers", () => {
        const { chart, wrapper } = fakeChart({
            spec: categorical,
            colorLegend: { display: true, pos: [10, 20], size: [0, "90"] },
        });
        chart.setColorLegend();
        expect(renderedLegend(wrapper).position.size).toBeUndefined();
    });

    test("keeps the size the user chose when the legend is redrawn", () => {
        const existingLegend = document.createElement("div");
        existingLegend.style.left = "30px";
        existingLegend.style.top = "40px";
        const { chart, wrapper } = fakeChart({ spec: categorical, existingLegend, userSize: [130, 70] });
        chart.setColorLegend();
        expect(renderedLegend(wrapper).position).toEqual({ left: "30px", top: "40px", size: [130, 70] });
    });

    test("saves the legend's size with its position", () => {
        const legend = document.createElement("div");
        legend.style.left = "448px";
        legend.style.top = "497px";
        const { chart } = fakeChart({ spec: categorical, existingLegend: legend, userSize: [142, 162] });
        expect(chart.getConfig().color_legend).toEqual({ display: true, pos: [448, 497], size: [142, 162] });
    });

    test("doesn't add a size the user never set", () => {
        const legend = document.createElement("div");
        const { chart } = fakeChart({ spec: categorical, existingLegend: legend });
        expect(chart.getConfig().color_legend).not.toHaveProperty("size");
    });

    test("saves the size when the legend is hidden", () => {
        const legend = document.createElement("div");
        legend.style.left = "5px";
        legend.style.top = "6px";
        const { chart } = fakeChart({
            spec: categorical,
            colorLegend: { display: false },
            existingLegend: legend,
            userSize: [100, 60],
        });
        chart.setColorLegend();
        expect(chart.config.color_legend).toEqual({ display: false, pos: [5, 6], size: [100, 60] });
    });
});

describe("BaseChart continuous color legend", () => {
    test("can be resized", () => {
        const { chart, wrapper } = fakeChart({ spec: continuous });
        chart.setColorLegend();
        expect(renderedLegend(wrapper).interactions.resizable).toBe(true);
    });

    test("has no range brush", () => {
        const { chart, wrapper } = fakeChart({ spec: continuous });
        chart.setColorLegend();
        const { legendProps } = renderedLegend(wrapper);
        expect(legendProps.onContinuousRangeChange).toBeUndefined();
        expect(legendProps.activeContinuousRange).toBeUndefined();
    });

    test("doesn't restore a range filter saved in an older view", () => {
        const { chart, getDimension } = fakeChart({
            spec: continuous,
            colorLegend: {
                display: true,
                filter: { kind: "continuous", column: "cd14", range: [2, 5] },
            },
        });
        chart.setColorLegend();
        expect(getDimension).not.toHaveBeenCalled();
        expect(chart.colorLegendFilter).toBeNull();
    });
});

describe("BaseChart categorical color legend", () => {
    test("still filters when a category is clicked", () => {
        const { chart, wrapper } = fakeChart({ spec: categorical });
        chart.setColorLegend();
        expect(renderedLegend(wrapper).legendProps.onCategoricalItemClick).toBeTypeOf("function");
    });
});
