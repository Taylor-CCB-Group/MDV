import type DataStore from "@/datastore/DataStore";
import type Dimension from "@/datastore/Dimension";
import type { ColorLegendSpec } from "@/react/legend/color_legend/types";

// Keep these operations as utilities while the state model is small.
// Convert this module into a controller if filter modes and lifecycle state grow.
export type ColorLegendCategoricalFilter = {
    kind: "categorical";
    column: string;
    value: string;
};

export type ColorLegendContinuousFilter = {
    kind: "continuous";
    column: string;
    range: [number, number];
};

export type ColorLegendFilter =
    | ColorLegendCategoricalFilter
    | ColorLegendContinuousFilter;

export type ColorLegendFilterChart = {
    config: {
        color_legend?: {
            display?: boolean;
            filter?: ColorLegendFilter;
        };
    };
    dataStore: DataStore;
    colorLegendFilterDimension: Dimension | null;
    colorLegendFilterDimensionKind: ColorLegendFilter["kind"] | null;
    colorLegendFilter: ColorLegendFilter | null;
    setColorLegend(): void;
    updateResetButtonVisibility(): void;
};

function getColorLegendFilterDimension(
    chart: ColorLegendFilterChart,
    filter: ColorLegendFilter,
): Dimension {
    if (
        chart.colorLegendFilterDimension &&
        chart.colorLegendFilterDimensionKind !== filter.kind
    ) {
        chart.colorLegendFilterDimension.destroy(false);
        chart.colorLegendFilterDimension = null;
        chart.colorLegendFilterDimensionKind = null;
    }
    if (!chart.colorLegendFilterDimension) {
        switch (filter.kind) {
            case "categorical":
                chart.colorLegendFilterDimension =
                    chart.dataStore.getDimension("category_dimension");
                break;
            case "continuous":
                chart.colorLegendFilterDimension =
                    chart.dataStore.getDimension("range_dimension");
                break;
        }
        chart.colorLegendFilterDimensionKind = filter.kind;
    }
    return chart.colorLegendFilterDimension;
}

export function getActiveCategoricalColorLegendValue(
    chart: ColorLegendFilterChart,
    spec: ColorLegendSpec,
): string | null {
    if (
        spec.kind !== "categorical" ||
        chart.colorLegendFilter?.kind !== "categorical" ||
        chart.colorLegendFilter.column !== spec.column
    ) {
        return null;
    }
    return chart.colorLegendFilter.value;
}

export function getActiveContinuousColorLegendRange(
    chart: ColorLegendFilterChart,
    spec: ColorLegendSpec,
): [number, number] | null {
    if (
        spec.kind !== "continuous" ||
        chart.colorLegendFilter?.kind !== "continuous" ||
        chart.colorLegendFilter.column !== spec.column
    ) {
        return null;
    }
    return chart.colorLegendFilter.range;
}

export function clearColorLegendFilter(
    chart: ColorLegendFilterChart,
    updateLegend = true,
): void {
    if (!chart.colorLegendFilter && !chart.config.color_legend?.filter) {
        return;
    }
    chart.colorLegendFilterDimension?.removeFilter();
    chart.colorLegendFilter = null;
    if (chart.config.color_legend) {
        delete chart.config.color_legend.filter;
    }
    chart.updateResetButtonVisibility();
    if (updateLegend) {
        chart.setColorLegend();
    }
}

export function applyColorLegendFilter(
    chart: ColorLegendFilterChart,
    filter: ColorLegendFilter,
    updateLegend = true,
): void {
    const dimension = getColorLegendFilterDimension(chart, filter);
    chart.colorLegendFilter = filter;
    chart.config.color_legend ??= { display: true };
    chart.config.color_legend.filter = filter;

    switch (filter.kind) {
        case "categorical":
            dimension.filter(
                "filterCategories",
                [filter.column],
                filter.value,
            );
            break;
        case "continuous":
            // This mirrors SelectionDialog's numeric filter path:
            // RangeDimension.filterRange keeps rows where column is between min and max.
            dimension.filter(
                "filterRange",
                [filter.column],
                { min: filter.range[0], max: filter.range[1] },
                true,
            );
            break;
    }
    chart.updateResetButtonVisibility();
    if (updateLegend) {
        chart.setColorLegend();
    }
}

export function toggleCategoricalColorLegendFilter(
    chart: ColorLegendFilterChart,
    spec: ColorLegendSpec,
    value: string,
): void {
    if (spec.kind !== "categorical") {
        return;
    }
    const isActive =
        chart.colorLegendFilter?.kind === "categorical" &&
        chart.colorLegendFilter.column === spec.column &&
        chart.colorLegendFilter.value === value;

    if (isActive) {
        clearColorLegendFilter(chart);
        return;
    }

    applyColorLegendFilter(chart, {
        kind: "categorical",
        column: spec.column,
        value,
    });
}

export function setContinuousColorLegendFilter(
    chart: ColorLegendFilterChart,
    spec: ColorLegendSpec,
    range: [number, number] | null,
): void {
    if (spec.kind !== "continuous") {
        return;
    }
    if (range === null) {
        clearColorLegendFilter(chart);
        return;
    }
    applyColorLegendFilter(chart, {
        kind: "continuous",
        column: spec.column,
        range,
    });
}

export function restoreColorLegendFilter(
    chart: ColorLegendFilterChart,
    spec: ColorLegendSpec,
): void {
    if (chart.colorLegendFilter) {
        chart.updateResetButtonVisibility();
        return;
    }
    const savedFilter = chart.config.color_legend?.filter;
    if (
        !savedFilter ||
        savedFilter.column !== spec.column ||
        savedFilter.kind !== spec.kind
    ) {
        return;
    }
    applyColorLegendFilter(chart, savedFilter, false);
}

export function destroyColorLegendFilter(
    chart: ColorLegendFilterChart,
    notify?: boolean,
): void {
    chart.colorLegendFilterDimension?.destroy(notify);
    chart.colorLegendFilterDimension = null;
    chart.colorLegendFilterDimensionKind = null;
    chart.colorLegendFilter = null;
}
