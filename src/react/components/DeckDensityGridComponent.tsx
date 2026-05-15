import DeckGL from "@deck.gl/react";
import { OrthographicView, type Layer, type OrthographicViewState } from "@deck.gl/core";
import type { Matrix4 } from "@math.gl/core";
import { action } from "mobx";
import { useCallback, useLayoutEffect, useMemo, useRef } from "react";
import { getFieldColor } from "../fieldColorManager";
import useGateLayers from "../hooks/useGateLayers";
import { useChartArrayMetrics } from "../hooks/useChartArrayMetrics";
import { useChartArrayVisibleIndices } from "../hooks/useChartArrayVisibleIndices";
import {
    useChartID,
    useChartSize,
    useConfig,
    useFieldSpecs,
    useFilteredIndices,
    useParamColumns,
} from "../hooks";
import { isLoadedNumericContourField, useFieldContour, type DualContourLegacyConfig } from "../contour_state";
import type { ScatterPlotConfig2D } from "../scatter_state";
import { useRange, useSpatialLayers } from "../spatial_context";
import ChartArrayLayout, { type ChartArrayLayoutHandle } from "./ChartArrayLayout";
import { useChart } from "../context";
import {
    getDensityGridViewId,
    getDensityGridViewStates,
    hasUsableOrthographicViewState,
    matchesDensityGridView,
} from "./densityGridUtils";

type DensityGridConfig = ScatterPlotConfig2D & DualContourLegacyConfig;
type CloneableDeckLayer = Layer & {
    clone: (props: Record<string, unknown>) => Layer;
};

function cloneDeckLayer(layer: CloneableDeckLayer, props: Record<string, unknown>): Layer {
    return layer.clone(props);
}

function getSerializableViewState(viewState: OrthographicViewState): OrthographicViewState {
    return {
        target: viewState.target,
        zoom: viewState.zoom,
        minZoom: viewState.minZoom,
        maxZoom: viewState.maxZoom,
    };
}

const BBOX_PERCENTILE_LOW = 0.01;
const BBOX_PERCENTILE_HIGH = 0.99;

function percentile(sorted: number[], p: number) {
    if (sorted.length === 0) return Number.NaN;
    const index = Math.min(sorted.length - 1, Math.max(0, Math.floor(p * (sorted.length - 1))));
    return sorted[index];
}

function fitViewStateToFilteredRows(
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

export default function DeckDensityGridComponent() {
    const chartId = useChartID();
    const chart = useChart() as { pendingRecenter?: boolean };
    const lastParamKeyRef = useRef("");
    const lastFilteredRowCountRef = useRef(-1);
    const scrollContainerRef = useRef<HTMLDivElement>(null);
    const layoutRef = useRef<ChartArrayLayoutHandle>(null);
    const config = useConfig<DensityGridConfig>();
    const [chartWidth, chartHeight] = useChartSize();
    const rows = useFilteredIndices();
    const {
        scatterProps: { scatterplotLayer, greyScatterplotLayer, setScatterKeyboardActive },
    } = useSpatialLayers();
    const [cx, cy] = useParamColumns();
    const { modelMatrix } = useRange();
    const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();
    const loadedFields = useFieldSpecs(config.densityFields);
    const densityFields = loadedFields.filter(isLoadedNumericContourField);
    const configuredFieldCount = Array.isArray(config.densityFields) ? config.densityFields.length : 0;
    const cellCount = densityFields.length;

    const metrics = useChartArrayMetrics(layoutRef, cellCount);
    const visibleCellIndices = useChartArrayVisibleIndices(scrollContainerRef, layoutRef, cellCount);
    const visibleFieldIndicesKey = visibleCellIndices.join("\u0000");
    const visibleFieldIndices = useMemo(
        () => visibleCellIndices,
        [visibleFieldIndicesKey],
    );

    const coordKey = `${cx.field}--${cy.field}`.replace(/[^A-Za-z0-9_-]/g, "_");
    const contourLayers = useFieldContour({
        id: `density-grid-${coordKey}`,
        fields: densityFields,
        visibleFieldIndices,
        fill: config.contour_fill,
        bandwidth: config.contour_bandwidth ?? 0.1,
        intensity: config.contour_intensity ?? 0.1,
        opacity: config.contour_opacity ?? 0.2,
        fillThreshold: config.contour_fillThreshold ?? 2,
    });
    const viewIds = useMemo(
        () => densityFields.map((field, index) => getDensityGridViewId(chartId, field.field, index)),
        [chartId, densityFields],
    );
    const visibleViewIds = useMemo(
        () => visibleCellIndices.map((index) => viewIds[index]).filter((viewId): viewId is string => Boolean(viewId)),
        [visibleCellIndices, viewIds],
    );
    const views = useMemo(
        () =>
            visibleCellIndices.map((index) => {
                const bounds = metrics.cellBounds[index];
                if (!bounds || bounds.width <= 0 || bounds.height <= 0) return null;
                return new OrthographicView({
                    id: viewIds[index],
                    x: bounds.x,
                    y: bounds.y,
                    width: bounds.width,
                    height: bounds.height,
                    flipY: false,
                    controller: { dragPan: controllerOptions.dragPan },
                });
            }).filter((view): view is OrthographicView => view !== null),
        [visibleCellIndices, metrics.cellBounds, viewIds, controllerOptions.dragPan],
    );
    const viewState = useMemo(
        () => getDensityGridViewStates(visibleViewIds, config.viewState),
        [visibleViewIds, config.viewState],
    );

    const paramKey = `${cx.field}\u0000${cy.field}`;
    const referenceCellWidth = metrics.cellBounds[0]?.width ?? 0;
    const referenceCellHeight = metrics.cellBounds[0]?.height ?? 0;

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

    const radiusPixels = useMemo(
        () => Math.max(1, 30 * (config.contour_bandwidth ?? 0.1) * 2 ** Number(config.viewState.zoom ?? 0)),
        [config.contour_bandwidth, config.viewState.zoom],
    );
    const layers = useMemo(
        () =>
            visibleCellIndices.flatMap((index) => {
                const contourLayer = contourLayers[index];
                const viewId = viewIds[index];
                if (!contourLayer || !viewId) return [];
                const greyLayer = cloneDeckLayer(greyScatterplotLayer, {
                    id: `${viewId}-grey`,
                    viewId,
                });
                const scatterLayer = cloneDeckLayer(scatterplotLayer, {
                    id: viewId,
                    viewId,
                    contourLayers: [
                        {
                            ...contourLayer,
                            radiusPixels,
                            debounce: 250,
                            weightsTextureSize: 256,
                        },
                    ],
                });
                const gateLayers = [
                    gateDisplayLayer
                        ? cloneDeckLayer(gateDisplayLayer, {
                              id: `${viewId}-gate`,
                              viewId,
                          })
                        : null,
                    gateLabelLayer
                        ? cloneDeckLayer(gateLabelLayer, {
                              id: `${viewId}-gate-label`,
                              viewId,
                          })
                        : null,
                ].filter((layer) => layer !== null);
                return [greyLayer, scatterLayer, ...gateLayers];
            }),
        [
            visibleCellIndices,
            contourLayers,
            viewIds,
            radiusPixels,
            scatterplotLayer,
            greyScatterplotLayer,
            gateDisplayLayer,
            gateLabelLayer,
            coordKey,
        ],
    );

    const renderCell = useCallback(
        (index: number) => {
            const field = densityFields[index];
            if (!field) return null;
            const color = getFieldColor(field.field);
            return (
                <div
                    className="pointer-events-none absolute left-2 top-2 z-10 max-w-[calc(100%-16px)] truncate rounded-sm px-2 py-1 text-xs"
                    style={{
                        backgroundColor: "rgba(15, 23, 42, 0.76)",
                        color: `rgb(${color[0]}, ${color[1]}, ${color[2]})`,
                    }}
                    title={field.name}
                >
                    {field.name || field.field}
                </div>
            );
        },
        [densityFields],
    );

    const hasValidViews =
        visibleCellIndices.length > 0 &&
        visibleCellIndices.every((index) => {
            const bounds = metrics.cellBounds[index];
            return (
                bounds !== undefined &&
                bounds.width > 0 &&
                bounds.height > 0 &&
                Number.isFinite(bounds.width) &&
                Number.isFinite(bounds.height)
            );
        });
    const hasCanvas =
        metrics.rootSize.width > 0 &&
        metrics.rootSize.height > 0 &&
        hasValidViews &&
        Number.isFinite(radiusPixels) &&
        radiusPixels > 0;

    const deckOverlay = useMemo(
        () =>
            hasCanvas ? (
                <DeckGL
                    controller={false}
                    layerFilter={({ layer, viewport }) => matchesDensityGridView(layer.id, viewport.id)}
                    layers={layers}
                    views={views}
                    viewState={viewState}
                    useDevicePixels={true}
                    onViewStateChange={({ viewState: nextViewState }: { viewState: OrthographicViewState }) => {
                        action(() => {
                            config.viewState = getSerializableViewState(nextViewState);
                        })();
                    }}
                    getCursor={({ isDragging }) => (isDragging ? "grabbing" : "crosshair")}
                />
            ) : null,
        [hasCanvas, layers, views, viewState, config],
    );

    if (configuredFieldCount === 0) {
        return <div className="flex h-full items-center justify-center text-sm">Choose density fields to build the grid.</div>;
    }
    if (densityFields.length === 0) {
        return <div className="flex h-full items-center justify-center text-sm">Loading density fields...</div>;
    }
    if (rows.length === 0) {
        return <div className="flex h-full items-center justify-center text-sm">No rows remain after the current filters.</div>;
    }

    return (
        <div className="relative flex h-full w-full min-w-0 flex-col">
            <div
                ref={scrollContainerRef}
                className="min-h-0 min-w-0 flex-1 overflow-auto overflow-x-hidden"
                onMouseDown={() => setScatterKeyboardActive(true)}
                onMouseEnter={() => setScatterKeyboardActive(true)}
                onMouseLeave={() => setScatterKeyboardActive(false)}
            >
                <ChartArrayLayout
                    cellCount={cellCount}
                    cellKeys={densityFields.map((field) => field.field)}
                    layoutRef={layoutRef}
                    canvasOverlay={deckOverlay}
                    renderCell={renderCell}
                />
            </div>
        </div>
    );
}
