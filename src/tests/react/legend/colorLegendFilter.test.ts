import { describe, expect, test, vi } from "vitest";
import DataStore from "@/datastore/DataStore";
import Dimension from "@/datastore/Dimension";
import {
    applyColorLegendFilter,
    clearColorLegendFilter,
    getActiveContinuousColorLegendRange,
    restoreColorLegendFilter,
    setContinuousColorLegendFilter,
    type ColorLegendFilterChart,
} from "@/react/legend/color_legend/colorLegendFilter";
import type { ColorLegendContinuousSpec } from "@/react/legend/color_legend/types";

function createChart(): ColorLegendFilterChart {
    const dimension = Object.assign(Object.create(Dimension.prototype), {
        filter: vi.fn(),
        removeFilter: vi.fn(),
        destroy: vi.fn(),
    });
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

const continuousSpec: ColorLegendContinuousSpec = {
    kind: "continuous",
    label: "Expression",
    column: "expression",
    colors: ["#000000", "#ffffff"],
    range: [0, 100],
};

describe("colorLegendFilter", () => {
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
