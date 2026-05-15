import DeckGL from "@deck.gl/react";
import { OrthographicView, type OrthographicViewState } from "@deck.gl/core";
import { useVirtualizer } from "@tanstack/react-virtual";
import { action } from "mobx";
import { useEffect, useMemo, useRef } from "react";
import { getFieldColor } from "../fieldColorManager";
import useGateLayers from "../hooks/useGateLayers";
import { useChartID, useChartSize, useConfig, useFieldSpecs, useFilteredIndices } from "../hooks";
import { isLoadedNumericContourField, useFieldContour, type DualContourLegacyConfig } from "../contour_state";
import type { ScatterPlotConfig2D } from "../scatter_state";
import { useSpatialLayers } from "../spatial_context";
import {
    getDensityGridCellBounds,
    getDensityGridLayout,
    getDensityGridViewId,
    getDensityGridViewStates,
    getDensityGridVisibleCellIndices,
    matchesDensityGridView,
} from "./densityGridUtils";

type DensityGridConfig = ScatterPlotConfig2D & DualContourLegacyConfig;
type CloneableDeckLayer = {
    clone: (props: any) => any;
};

function cloneDeckLayer(layer: CloneableDeckLayer, props: Record<string, unknown>) {
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

const DENSITY_GRID_VIRTUAL_OVERSCAN = 1;

export default function DeckDensityGridComponent() {
    const chartId = useChartID();
    const initializedViewRef = useRef(false);
    const scrollContainerRef = useRef<HTMLDivElement>(null);
    const [viewportWidth, viewportHeight] = useChartSize();
    const config = useConfig<DensityGridConfig>();
    const rows = useFilteredIndices();
    const {
        scatterProps: { scatterplotLayer, greyScatterplotLayer, setScatterKeyboardActive },
    } = useSpatialLayers();
    const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();
    const loadedFields = useFieldSpecs(config.densityFields);
    const densityFields = loadedFields.filter(isLoadedNumericContourField);
    const configuredFieldCount = Array.isArray(config.densityFields) ? config.densityFields.length : 0;

    const layout = useMemo(
        () => getDensityGridLayout(viewportWidth, viewportHeight, Math.max(densityFields.length, 1)),
        [viewportWidth, viewportHeight, densityFields.length],
    );
    const rowVirtualizer = useVirtualizer({
        count: layout.rows,
        getScrollElement: () => scrollContainerRef.current,
        estimateSize: () => layout.cellHeight,
        overscan: DENSITY_GRID_VIRTUAL_OVERSCAN,
    });
    const columnVirtualizer = useVirtualizer({
        count: layout.columns,
        getScrollElement: () => scrollContainerRef.current,
        estimateSize: () => layout.cellWidth,
        horizontal: true,
        overscan: DENSITY_GRID_VIRTUAL_OVERSCAN,
    });
    const visibleCellIndices = getDensityGridVisibleCellIndices(
        layout,
        densityFields.length,
        rowVirtualizer.getVirtualItems(),
        columnVirtualizer.getVirtualItems(),
    );
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
                const bounds = getDensityGridCellBounds(layout, index);
                return new OrthographicView({
                    id: viewIds[index],
                    x: bounds.x,
                    y: bounds.y,
                    width: bounds.width,
                    height: bounds.height,
                    flipY: false,
                    controller: { dragPan: controllerOptions.dragPan },
                });
            }),
        [visibleCellIndices, layout, viewIds, controllerOptions.dragPan],
    );
    const viewState = useMemo(
        () => getDensityGridViewStates(visibleViewIds, config.viewState),
        [visibleViewIds, config.viewState],
    );
    useEffect(() => {
        if (initializedViewRef.current) return;
        if (densityFields.length === 0) return;
        initializedViewRef.current = true;
        action(() => {
            config.viewState = {
                ...config.viewState,
                zoom:
                    Number(config.viewState.zoom ?? 0) +
                    Math.log2(
                        Math.max(
                            Number.MIN_VALUE,
                            layout.cellSize / Math.max(viewportWidth, 1),
                        ),
                    ),
            };
        })();
    }, [config, densityFields.length, layout.cellSize, viewportWidth]);
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
        <div className="relative flex h-full w-full flex-col">
            <div
                ref={scrollContainerRef}
                className="min-h-0 flex-1 overflow-auto"
                onMouseDown={() => setScatterKeyboardActive(true)}
                onMouseEnter={() => setScatterKeyboardActive(true)}
                onMouseLeave={() => setScatterKeyboardActive(false)}
            >
                <div
                    className="relative"
                    style={{ width: layout.contentWidth, height: layout.contentHeight }}
                >
                    <DeckGL
                        controller={false}
                        width={layout.contentWidth}
                        height={layout.contentHeight}
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
                    <div
                        className="pointer-events-none absolute inset-0"
                        style={{ width: layout.contentWidth, height: layout.contentHeight }}
                    >
                        {visibleCellIndices.map((index) => {
                            const field = densityFields[index];
                            if (!field) return null;
                            const bounds = getDensityGridCellBounds(layout, index);
                            const color = getFieldColor(field.field);
                            return (
                                <div
                                    key={field.field}
                                    className="absolute border border-slate-400/25"
                                    style={{
                                        left: `${bounds.x}px`,
                                        top: `${bounds.y}px`,
                                        width: `${bounds.width}px`,
                                        height: `${bounds.height}px`,
                                    }}
                                >
                                    <div
                                        className="absolute left-2 top-2 max-w-[calc(100%-16px)] truncate rounded-sm px-2 py-1 text-xs"
                                        style={{
                                            backgroundColor: "rgba(15, 23, 42, 0.76)",
                                            color: `rgb(${color[0]}, ${color[1]}, ${color[2]})`,
                                        }}
                                        title={field.name}
                                    >
                                        {field.name || field.field}
                                    </div>
                                </div>
                            );
                        })}
                    </div>
                </div>
            </div>
        </div>
    );
}
