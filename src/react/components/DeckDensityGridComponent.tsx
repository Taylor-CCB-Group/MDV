import DeckGL from "@deck.gl/react";
import { OrthographicView, type OrthographicViewState } from "@deck.gl/core";
import { action } from "mobx";
import { useEffect, useMemo, useRef } from "react";
import { getFieldColor } from "../fieldColorManager";
import { useChartID, useChartSize, useConfig, useFieldSpecs, useFilteredIndices } from "../hooks";
import { isLoadedNumericContourField, useFieldContour, type DualContourLegacyConfig } from "../contour_state";
import type { ScatterPlotConfig2D } from "../scatter_state";
import { useSpatialLayers } from "../spatial_context";
import {
    getDensityGridCellBounds,
    getDensityGridLayout,
    getDensityGridViewId,
    getDensityGridViewStates,
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

export default function DeckDensityGridComponent() {
    const chartId = useChartID();
    const initializedViewRef = useRef(false);
    const [width, height] = useChartSize();
    const config = useConfig<DensityGridConfig>();
    const rows = useFilteredIndices();
    const {
        scatterProps: { scatterplotLayer, greyScatterplotLayer },
    } = useSpatialLayers();
    const loadedFields = useFieldSpecs(config.densityFields);
    const densityFields = loadedFields.filter(isLoadedNumericContourField);
    const configuredFieldCount = Array.isArray(config.densityFields) ? config.densityFields.length : 0;
    const contourLayers = useFieldContour({
        id: "density-grid",
        fields: densityFields,
        fill: config.contour_fill,
        bandwidth: config.contour_bandwidth ?? 0.1,
        intensity: config.contour_intensity ?? 0.1,
        opacity: config.contour_opacity ?? 0.2,
        fillThreshold: config.contour_fillThreshold ?? 2,
    });

    const layout = useMemo(
        () => getDensityGridLayout(width, height, Math.max(densityFields.length, 1)),
        [width, height, densityFields.length],
    );
    const viewIds = useMemo(
        () => densityFields.map((field, index) => getDensityGridViewId(chartId, field.field, index)),
        [chartId, densityFields],
    );
    const views = useMemo(
        () =>
            densityFields.map((field, index) => {
                const bounds = getDensityGridCellBounds(layout, index);
                return new OrthographicView({
                    id: viewIds[index],
                    x: bounds.x,
                    y: bounds.y,
                    width: bounds.width,
                    height: bounds.height,
                    flipY: false,
                    // Top-level DeckGL `controller` only applies to the first view.
                    controller: { dragPan: true },
                });
            }),
        [densityFields, layout, viewIds],
    );
    const viewState = useMemo(
        () => getDensityGridViewStates(viewIds, config.viewState),
        [viewIds, config.viewState],
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
                            Math.min(layout.cellWidth / Math.max(width, 1), layout.cellHeight / Math.max(height, 1)),
                        ),
                    ),
            };
        })();
    }, [config, densityFields.length, height, layout.cellHeight, layout.cellWidth, width]);
    const radiusPixels = useMemo(
        () => Math.max(1, 30 * (config.contour_bandwidth ?? 0.1) * 2 ** Number(config.viewState.zoom ?? 0)),
        [config.contour_bandwidth, config.viewState.zoom],
    );
    const layers = useMemo(
        () =>
            contourLayers.flatMap((contourLayer, index) => {
                const viewId = viewIds[index];
                if (!viewId) return [];
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
                return [greyLayer, scatterLayer];
            }),
        [contourLayers, viewIds, radiusPixels, scatterplotLayer, greyScatterplotLayer],
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
        <div className="relative h-full w-full">
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
            <div className="pointer-events-none absolute inset-0">
                {densityFields.map((field, index) => {
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
    );
}
