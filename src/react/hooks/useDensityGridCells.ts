import type { Matrix4 } from "@math.gl/core";
import type { OrthographicViewState } from "@deck.gl/core";
import { action } from "mobx";
import { useLayoutEffect, useMemo, useRef } from "react";
import {
    isLoadedNumericContourField,
    useFieldContour,
    type DualContourLegacyConfig,
} from "../contour_state";
import { useChart } from "../context";
import { useChartSize, useConfig, useFieldSpecs, useFilteredIndices, useParamColumns } from "../hooks";
import type { LoadedDataColumn, NumberDataType } from "@/charts/charts";
import { useRange } from "../spatial_context";
import type { ScatterPlotConfig2D } from "../scatter_state";
import {
    hasUsableOrthographicViewState,
} from "../components/chartArrayGridUtils";
import type { ChartArrayCell } from "./useChartArrayGrid";

type DensityGridConfig = ScatterPlotConfig2D & DualContourLegacyConfig;

const BBOX_PERCENTILE_LOW = 0.01;
const BBOX_PERCENTILE_HIGH = 0.99;

function percentile(sorted: number[], p: number) {
    if (sorted.length === 0) return Number.NaN;
    const index = Math.min(sorted.length - 1, Math.max(0, Math.floor(p * (sorted.length - 1))));
    return sorted[index];
}

export function fitViewStateToFilteredRows(
    rows: Uint32Array,
    cx: { data: ArrayLike<number> },
    cy: { data: ArrayLike<number> },
    modelMatrix: Matrix4,
    viewportWidth: number,
    viewportHeight: number,
): Pick<OrthographicViewState, "target" | "zoom"> | null {
    if (rows.length === 0 || viewportWidth <= 0 || viewportHeight <= 0) return null;

    const xs: number[] = [];
    const ys: number[] = [];
    for (let i = 0; i < rows.length; i++) {
        const row = rows[i];
        const x = cx.data[row];
        const y = cy.data[row];
        if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
        xs.push(x);
        ys.push(y);
    }
    if (xs.length === 0 || ys.length === 0) return null;

    xs.sort((a, b) => a - b);
    ys.sort((a, b) => a - b);
    let minX = percentile(xs, BBOX_PERCENTILE_LOW);
    let maxX = percentile(xs, BBOX_PERCENTILE_HIGH);
    let minY = percentile(ys, BBOX_PERCENTILE_LOW);
    let maxY = percentile(ys, BBOX_PERCENTILE_HIGH);
    if (!Number.isFinite(minX) || !Number.isFinite(maxX) || !Number.isFinite(minY) || !Number.isFinite(maxY)) {
        return null;
    }

    const minWorld = modelMatrix.transformAsPoint([minX, minY, 0]);
    const maxWorld = modelMatrix.transformAsPoint([maxX, maxY, 0]);
    minX = minWorld[0];
    minY = minWorld[1];
    maxX = maxWorld[0];
    maxY = maxWorld[1];

    const dx = maxX - minX;
    const dy = maxY - minY;
    const epsilon = Math.max(1e-9, Math.min(Math.abs(minX), Math.abs(maxX)) * 1e-6);
    const safeDx = Math.max(epsilon, dx);
    const safeDy = Math.max(epsilon, dy);
    let zoom = Math.log2(Math.min(viewportWidth / safeDx, viewportHeight / safeDy)) - 0.6;
    if (!Number.isFinite(zoom)) zoom = 0;

    return {
        target: [(minX + maxX) / 2, (minY + maxY) / 2, 0],
        zoom,
    };
}

export function useDensityGridCells() {
    const config = useConfig<DensityGridConfig>();
    const loadedFields = useFieldSpecs(config.densityFields);
    const densityFields = loadedFields.filter(isLoadedNumericContourField);
    const configuredFieldCount = Array.isArray(config.densityFields) ? config.densityFields.length : 0;

    const cells: ChartArrayCell[] = useMemo(
        () =>
            densityFields.map((field) => ({
                key: field.field,
                label: field.name || field.field,
            })),
        [densityFields],
    );

    return {
        cells,
        densityFields,
        configuredFieldCount,
    };
}

export function useDensityGridContours(
    densityFields: LoadedDataColumn<NumberDataType>[],
    visibleFieldIndices: readonly number[],
) {
    const config = useConfig<DensityGridConfig>();
    const [cx, cy] = useParamColumns();
    const coordKey = `${cx.field}--${cy.field}`.replace(/[^A-Za-z0-9_-]/g, "_");

    return useFieldContour({
        id: `density-grid-${coordKey}`,
        fields: densityFields,
        visibleFieldIndices,
        fill: config.contour_fill,
        bandwidth: config.contour_bandwidth ?? 0.1,
        intensity: config.contour_intensity ?? 0.1,
        opacity: config.contour_opacity ?? 0.2,
        fillThreshold: config.contour_fillThreshold ?? 2,
    });
}

export function useDensityGridViewState(
    referenceCellWidth: number,
    referenceCellHeight: number,
) {
    const chart = useChart() as { pendingRecenter?: boolean };
    const config = useConfig<DensityGridConfig>();
    const [chartWidth, chartHeight] = useChartSize();
    const rows = useFilteredIndices();
    const [cx, cy] = useParamColumns();
    const { modelMatrix } = useRange();
    const lastParamKeyRef = useRef("");
    const lastFilteredRowCountRef = useRef(-1);
    const paramKey = `${cx.field}\u0000${cy.field}`;

    useLayoutEffect(() => {
        const paramChanged = lastParamKeyRef.current !== "" && lastParamKeyRef.current !== paramKey;
        lastParamKeyRef.current = paramKey;

        const filterChanged = lastFilteredRowCountRef.current !== rows.length;
        lastFilteredRowCountRef.current = rows.length;

        const viewportWidth = referenceCellWidth > 0 ? referenceCellWidth : chartWidth;
        const viewportHeight = referenceCellHeight > 0 ? referenceCellHeight : chartHeight;
        const pendingRecenter = !!chart.pendingRecenter;
        const needsFullFit =
            pendingRecenter ||
            paramChanged ||
            !hasUsableOrthographicViewState(config.viewState) ||
            (config.zoom_on_filter && filterChanged);

        if (needsFullFit && viewportWidth > 0 && viewportHeight > 0) {
            const fitted = fitViewStateToFilteredRows(
                rows,
                cx,
                cy,
                modelMatrix,
                viewportWidth,
                viewportHeight,
            );
            if (!fitted) return;
            action(() => {
                chart.pendingRecenter = false;
                config.viewState = {
                    ...config.viewState,
                    target: fitted.target,
                    zoom: fitted.zoom,
                };
            })();
        }
    }, [
        paramKey,
        rows.length,
        cx,
        cy,
        modelMatrix,
        referenceCellWidth,
        referenceCellHeight,
        chartWidth,
        chartHeight,
        config,
        chart,
    ]);

    return useMemo(
        () => Math.max(1, 30 * (config.contour_bandwidth ?? 0.1) * 2 ** Number(config.viewState.zoom ?? 0)),
        [config.contour_bandwidth, config.viewState.zoom],
    );
}
