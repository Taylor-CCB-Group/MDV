import type { OrthographicViewState } from "@deck.gl/core";
import type { CategoricalDataType, LoadedDataColumn, NumberDataType } from "../../charts/charts";
import type { BaseConfig } from "../../charts/BaseChart";
import type { FieldSpec, FieldSpecs } from "@/lib/columnTypeHelpers";
import type { ContourVisualConfig } from "../contour_state";

export const splatterDefaults = {
    contour_fill: true,
    contour_bandwidth: 0.1,
    contour_intensity: 1,
    contour_opacity: 0.35,
    contour_fillThreshold: 2,
} satisfies ContourVisualConfig;

export type SplatterPlotConfig = BaseConfig &
    ContourVisualConfig & {
        category?: FieldSpec;
        densityFields?: FieldSpecs;
    };

export type VisibleSplatterCategory = {
    categoryIndex: number;
    label: string;
    rows: Uint32Array;
};

export type SplatterLayout = {
    headerHeight: number;
    labelWidth: number;
    plotWidth: number;
    plotHeight: number;
    cellWidth: number;
    cellHeight: number;
};

export function adaptSplatterConfig(originalConfig: SplatterPlotConfig): SplatterPlotConfig {
    const config = {
        ...splatterDefaults,
        ...originalConfig,
    };
    if (Array.isArray(config.param) && config.param.length > 2 && !config.densityFields) {
        config.densityFields = config.param.slice(2);
        config.param = config.param.slice(0, 2);
    }
    return config;
}

export function getVisibleSplatterCategories(
    categoryColumn: LoadedDataColumn<CategoricalDataType> | undefined,
    filteredRows: Uint32Array,
): VisibleSplatterCategory[] {
    if (!categoryColumn?.data) return [];
    const buckets = categoryColumn.values.map(() => [] as number[]);
    for (const rowIndex of filteredRows) {
        const categoryIndex = categoryColumn.data[rowIndex];
        if (categoryIndex === undefined) continue;
        buckets[categoryIndex]?.push(rowIndex);
    }
    const categories: VisibleSplatterCategory[] = [];
    for (const [categoryIndex, rows] of buckets.entries()) {
        if (rows.length === 0) continue;
        categories.push({
            categoryIndex,
            label: categoryColumn.values[categoryIndex] ?? "",
            rows: Uint32Array.from(rows),
        });
    }
    return categories;
}

export function getSplatterLayout(
    width: number,
    height: number,
    categories: string[],
    fields: string[],
): SplatterLayout {
    const longestCategory = categories.reduce((max, label) => Math.max(max, label.length), 0);
    const longestField = fields.reduce((max, label) => Math.max(max, label.length), 0);
    const labelWidth = Math.min(240, Math.max(88, longestCategory * 7 + 24));
    const headerHeight = Math.min(110, Math.max(58, longestField * 4 + 28));
    const plotWidth = Math.max(1, width - labelWidth - 8);
    const plotHeight = Math.max(1, height - headerHeight - 8);
    return {
        headerHeight,
        labelWidth,
        plotWidth,
        plotHeight,
        cellWidth: plotWidth / Math.max(fields.length, 1),
        cellHeight: plotHeight / Math.max(categories.length, 1),
    };
}

export function getSharedSplatterViewState(
    cx: LoadedDataColumn<NumberDataType>,
    cy: LoadedDataColumn<NumberDataType>,
    rows: Uint32Array,
    cellWidth: number,
    cellHeight: number,
): OrthographicViewState {
    if (rows.length === 0) {
        return {
            target: [0, 0, 0],
            zoom: 0,
            minZoom: -50,
        };
    }

    let minX = Number.POSITIVE_INFINITY;
    let maxX = Number.NEGATIVE_INFINITY;
    let minY = Number.POSITIVE_INFINITY;
    let maxY = Number.NEGATIVE_INFINITY;

    for (const rowIndex of rows) {
        const x = cx.data[rowIndex];
        const y = cy.data[rowIndex];
        if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
    }

    if (!Number.isFinite(minX) || !Number.isFinite(minY) || !Number.isFinite(maxX) || !Number.isFinite(maxY)) {
        return {
            target: [0, 0, 0],
            zoom: 0,
            minZoom: -50,
        };
    }

    const dx = maxX - minX;
    const dy = maxY - minY;
    const epsilon =
        Math.max(1e-9, Math.min(Math.abs(minX), Math.abs(maxX), Math.abs(minY), Math.abs(maxY), 1) * 1e-6);
    const safeDx = Math.max(dx, epsilon);
    const safeDy = Math.max(dy, epsilon);
    const paddedWidth = Math.max(cellWidth * 0.82, 1);
    const paddedHeight = Math.max(cellHeight * 0.82, 1);

    return {
        target: [(minX + maxX) / 2, (minY + maxY) / 2, 0],
        zoom: Math.log2(Math.min(paddedWidth / safeDx, paddedHeight / safeDy)),
        minZoom: -50,
    };
}
