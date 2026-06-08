import { describe, expect, test, vi } from "vitest";
import DataStore from "@/datastore/DataStore";
import Dimension from "@/datastore/Dimension";
import {
    applyColorLegendFilter,
    clearColorLegendFilter,
    getActiveCategoricalColorLegendValue,
    getActiveContinuousColorLegendRange,
    restoreColorLegendFilter,
    setContinuousColorLegendFilter,
    toggleCategoricalColorLegendFilter,
    type ColorLegendFilterChart,
} from "@/react/legend/color_legend/colorLegendFilter";
import type {
    ColorLegendCategoricalSpec,
    ColorLegendContinuousSpec,
} from "@/react/legend/color_legend/types";

function createChart(): ColorLegendFilterChart {
    const dimension = createDimension();
    const dataStore = Object.assign(Object.create(DataStore.prototype), {
        getDimension: vi.fn(() => dimension),
    });

    return {
        config: {},
        dataStore,
        colorLegendFilterDimension: null,
        colorLegendFilterDimensionKind: null,
        colorLegendFilter: null,
        setColorLegend: vi.fn(),
        updateResetButtonVisibility: vi.fn(),
    };
}

function createDimension() {
    return Object.assign(Object.create(Dimension.prototype), {
        filter: vi.fn(),
        removeFilter: vi.fn(),
        destroy: vi.fn(),
    });
}

function createChartWithDimensions(): ColorLegendFilterChart {
    const categoryDimension = createDimension();
    const rangeDimension = createDimension();
    const dataStore = Object.assign(Object.create(DataStore.prototype), {
        getDimension: vi.fn((type: string) =>
            type === "category_dimension" ? categoryDimension : rangeDimension,
        ),
    });

    return {
        config: {},
        dataStore,
        colorLegendFilterDimension: null,
        colorLegendFilterDimensionKind: null,
        colorLegendFilter: null,
        setColorLegend: vi.fn(),
        updateResetButtonVisibility: vi.fn(),
    };
}

const categoricalSpec: ColorLegendCategoricalSpec = {
    kind: "categorical",
    label: "Cell type",
    column: "cell_type",
    items: [
        {
            color: "#ff0000",
            name: "T-cell",
            value: "T-cell",
        },
        {
            color: "#00ff00",
            name: "B-cell",
            value: "B-cell",
        },
    ],
};

const continuousSpec: ColorLegendContinuousSpec = {
    kind: "continuous",
    label: "Expression",
    column: "expression",
    colors: ["#000000", "#ffffff"],
    range: [0, 100],
};

describe("colorLegendFilter", () => {
    test("applies categorical filters through category_dimension filterCategories", () => {
        const chart = createChart();

        applyColorLegendFilter(chart, {
            kind: "categorical",
            column: "cell_type",
            value: "T-cell",
        });

        expect(chart.dataStore.getDimension).toHaveBeenCalledWith("category_dimension");
        expect(chart.colorLegendFilterDimension?.filter).toHaveBeenCalledWith(
            "filterCategories",
            ["cell_type"],
            "T-cell",
        );
        expect(chart.colorLegendFilter).toEqual({
            kind: "categorical",
            column: "cell_type",
            value: "T-cell",
        });
        expect(chart.config.color_legend?.filter).toEqual(chart.colorLegendFilter);
        expect(chart.updateResetButtonVisibility).toHaveBeenCalled();
        expect(chart.setColorLegend).toHaveBeenCalled();
    });

    test("toggles categorical filters and exposes the active value", () => {
        const chart = createChart();

        toggleCategoricalColorLegendFilter(chart, categoricalSpec, "T-cell");
        expect(getActiveCategoricalColorLegendValue(chart, categoricalSpec)).toBe(
            "T-cell",
        );

        toggleCategoricalColorLegendFilter(chart, categoricalSpec, "T-cell");
        expect(chart.colorLegendFilter).toBeNull();
        expect(chart.config.color_legend?.filter).toBeUndefined();
        expect(chart.colorLegendFilterDimension?.removeFilter).toHaveBeenCalled();
    });

    test("replaces incompatible filter dimensions when switching filter kind", () => {
        const chart = createChartWithDimensions();

        applyColorLegendFilter(chart, {
            kind: "categorical",
            column: "cell_type",
            value: "T-cell",
        });
        const categoryDimension = chart.colorLegendFilterDimension;

        applyColorLegendFilter(chart, {
            kind: "continuous",
            column: "expression",
            range: [20, 60],
        });

        expect(categoryDimension?.destroy).toHaveBeenCalledWith(false);
        expect(chart.dataStore.getDimension).toHaveBeenCalledWith("range_dimension");
        expect(chart.colorLegendFilterDimensionKind).toBe("continuous");
        expect(chart.colorLegendFilterDimension?.filter).toHaveBeenCalledWith(
            "filterRange",
            ["expression"],
            { min: 20, max: 60 },
            true,
        );
    });

    test("applies continuous filters through range_dimension filterRange", () => {
        const chart = createChart();

        applyColorLegendFilter(chart, {
            kind: "continuous",
            column: "expression",
            range: [20, 60],
        });

        expect(chart.dataStore.getDimension).toHaveBeenCalledWith("range_dimension");
        expect(chart.colorLegendFilterDimension?.filter).toHaveBeenCalledWith(
            "filterRange",
            ["expression"],
            { min: 20, max: 60 },
            true,
        );
        expect(chart.colorLegendFilter).toEqual({
            kind: "continuous",
            column: "expression",
            range: [20, 60],
        });
        expect(chart.config.color_legend?.filter).toEqual(chart.colorLegendFilter);
        expect(chart.updateResetButtonVisibility).toHaveBeenCalled();
        expect(chart.setColorLegend).toHaveBeenCalled();
    });

    test("sets and clears continuous filter ranges", () => {
        const chart = createChart();

        setContinuousColorLegendFilter(chart, continuousSpec, [10, 50]);
        expect(getActiveContinuousColorLegendRange(chart, continuousSpec)).toEqual([
            10,
            50,
        ]);

        setContinuousColorLegendFilter(chart, continuousSpec, null);
        expect(chart.colorLegendFilter).toBeNull();
        expect(chart.config.color_legend?.filter).toBeUndefined();
        expect(chart.colorLegendFilterDimension?.removeFilter).toHaveBeenCalled();
    });

    test("restores saved continuous filters that match the current legend spec", () => {
        const chart = createChart();
        chart.config.color_legend = {
            display: true,
            filter: {
                kind: "continuous",
                column: "expression",
                range: [25, 75],
            },
        };

        restoreColorLegendFilter(chart, continuousSpec);

        expect(chart.colorLegendFilter).toEqual({
            kind: "continuous",
            column: "expression",
            range: [25, 75],
        });
        expect(chart.colorLegendFilterDimension?.filter).toHaveBeenCalledWith(
            "filterRange",
            ["expression"],
            { min: 25, max: 75 },
            true,
        );
        expect(chart.setColorLegend).not.toHaveBeenCalled();
    });

    test("clears continuous filters with the shared clear helper", () => {
        const chart = createChart();
        applyColorLegendFilter(chart, {
            kind: "continuous",
            column: "expression",
            range: [20, 60],
        });

        clearColorLegendFilter(chart);

        expect(chart.colorLegendFilter).toBeNull();
        expect(chart.config.color_legend?.filter).toBeUndefined();
        expect(chart.colorLegendFilterDimension?.removeFilter).toHaveBeenCalled();
    });
});
