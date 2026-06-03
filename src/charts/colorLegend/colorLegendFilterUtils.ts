import type DataStore from "@/datastore/DataStore";
import type Dimension from "@/datastore/Dimension";
import type { ColorLegendSpec } from "@/react/colorLegend/types";

// Keep these operations as utilities while the state model is small.
// Convert this module into a controller if filter modes and lifecycle state grow.
export type ColorLegendCategoricalFilter = {
    kind: "categorical";
    column: string;
    value: string;
};

// Add a numeric range variant to this union when continuous legend filtering is implemented.
export type ColorLegendFilter = ColorLegendCategoricalFilter;

export type ColorLegendFilterChart = {
    config: {
        color_legend?: {
            display?: boolean;
            filter?: ColorLegendFilter;
        };
    };
    dataStore: DataStore;
    colorLegendFilterDimension: Dimension | null;
    colorLegendFilter: ColorLegendFilter | null;
    setColorLegend(): void;
    updateResetButtonVisibility(): void;
};

function getColorLegendFilterDimension(
    chart: ColorLegendFilterChart,
    filter: ColorLegendFilter,
): Dimension {
    if (!chart.colorLegendFilterDimension) {
        switch (filter.kind) {
            case "categorical":
                chart.colorLegendFilterDimension =
                    chart.dataStore.getDimension("category_dimension");
                break;
            // Add a range dimension case for a future numeric legend filter variant.
        }
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
        // Add numeric range filtering here when continuous legend interactions are defined.
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
        savedFilter.kind !== "categorical" ||
        spec.kind !== "categorical" ||
        savedFilter.column !== spec.column
    ) {
        return;
    }
    // Restore numeric range filters here when continuous legend filtering is implemented.
    applyColorLegendFilter(chart, savedFilter, false);
}

export function destroyColorLegendFilter(
    chart: ColorLegendFilterChart,
    notify?: boolean,
): void {
    chart.colorLegendFilterDimension?.destroy(notify);
    chart.colorLegendFilterDimension = null;
    chart.colorLegendFilter = null;
}
