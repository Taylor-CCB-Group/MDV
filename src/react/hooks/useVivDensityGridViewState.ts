import { action } from "mobx";
import { useLayoutEffect, useMemo, useRef } from "react";
import { useChart } from "../context";
import { useChartSize, useConfig, useFilteredIndices, useParamColumns } from "../hooks";
import type { DualContourLegacyConfig } from "../contour_state";
import type { ScatterPlotConfig2D } from "../scatter_state";
import { useRange } from "../spatial_context";
import { useViewerStore, useViewerStoreApi } from "../components/avivatorish/state";
import { fitViewStateToFilteredRows } from "./useDensityGridCells";
import { hasUsableOrthographicViewState } from "../components/chartArrayGridUtils";

type VivDensityGridConfig = ScatterPlotConfig2D & DualContourLegacyConfig;

export function useVivDensityGridViewState(
    referenceCellWidth: number,
    referenceCellHeight: number,
) {
    const chart = useChart() as { pendingRecenter?: boolean };
    const config = useConfig<VivDensityGridConfig>();
    const [chartWidth, chartHeight] = useChartSize();
    const rows = useFilteredIndices();
    const [cx, cy] = useParamColumns();
    const { modelMatrix } = useRange();
    const viewerStore = useViewerStoreApi();
    const viewState = useViewerStore((store) => store.viewState);
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
            !hasUsableOrthographicViewState(viewState) ||
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
                viewerStore.setState({
                    viewState: {
                        ...viewState,
                        target: fitted.target,
                        zoom: fitted.zoom,
                    },
                });
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
        viewState,
        viewerStore,
    ]);

    return useMemo(
        () => Math.max(1, 30 * (config.contour_bandwidth ?? 0.1) * 2 ** Number(viewState?.zoom ?? 0)),
        [config.contour_bandwidth, viewState?.zoom],
    );
}
