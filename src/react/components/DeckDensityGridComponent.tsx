import DeckGL from "@deck.gl/react";
import { OrthographicView, type Layer, type OrthographicViewState } from "@deck.gl/core";
import { action } from "mobx";
import { useCallback, useEffect, useMemo, useRef } from "react";
import { getFieldColor } from "../fieldColorManager";
import useGateLayers from "../hooks/useGateLayers";
import { useChartArrayMetrics } from "../hooks/useChartArrayMetrics";
import { useChartArrayVisibleIndices } from "../hooks/useChartArrayVisibleIndices";
import { useChartID, useConfig, useFieldSpecs, useFilteredIndices } from "../hooks";
import { isLoadedNumericContourField, useFieldContour, type DualContourLegacyConfig } from "../contour_state";
import type { ScatterPlotConfig2D } from "../scatter_state";
import { useSpatialLayers } from "../spatial_context";
import ChartArrayLayout, { type ChartArrayLayoutHandle } from "./ChartArrayLayout";
import {
    getDensityGridViewId,
    getDensityGridViewStates,
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

export default function DeckDensityGridComponent() {
    const chartId = useChartID();
    const initializedViewRef = useRef(false);
    const scrollContainerRef = useRef<HTMLDivElement>(null);
    const layoutRef = useRef<ChartArrayLayoutHandle>(null);
    const config = useConfig<DensityGridConfig>();
    const rows = useFilteredIndices();
    const {
        scatterProps: { scatterplotLayer, greyScatterplotLayer, setScatterKeyboardActive },
    } = useSpatialLayers();
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

    const contourLayers = useFieldContour({
        id: "density-grid",
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

    const referenceCellWidth = metrics.cellBounds[0]?.width ?? 0;
    useEffect(() => {
        if (initializedViewRef.current) return;
        if (densityFields.length === 0 || referenceCellWidth <= 0) return;
        initializedViewRef.current = true;
        action(() => {
            config.viewState = {
                ...config.viewState,
                zoom:
                    Number(config.viewState.zoom ?? 0) +
                    Math.log2(
                        Math.max(
                            Number.MIN_VALUE,
                            referenceCellWidth / Math.max(metrics.rootSize.width, 1),
                        ),
                    ),
            };
        })();
    }, [config, densityFields.length, referenceCellWidth, metrics.rootSize.width]);

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
                            weightsTextureSize: 128,
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

    const { width: canvasWidth, height: canvasHeight } = metrics.rootSize;
    const hasCanvas = canvasWidth > 0 && canvasHeight > 0;

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
                <div className="relative">
                    {hasCanvas && (
                        <div
                            className="pointer-events-auto absolute left-0 top-0"
                            style={{ width: canvasWidth, height: canvasHeight }}
                        >
                            <DeckGL
                                controller={false}
                                width={canvasWidth}
                                height={canvasHeight}
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
                        </div>
                    )}
                    <ChartArrayLayout
                        cellCount={cellCount}
                        cellKeys={densityFields.map((field) => field.field)}
                        layoutRef={layoutRef}
                        className="pointer-events-none relative z-10"
                        renderCell={renderCell}
                    />
                </div>
            </div>
        </div>
    );
}
