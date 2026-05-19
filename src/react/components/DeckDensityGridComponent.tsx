import DeckGL from "@deck.gl/react";
import { OrthographicView, type OrthographicViewState } from "@deck.gl/core";
import { action } from "mobx";
import { useCallback, useMemo, useRef } from "react";
import { getFieldColor } from "../fieldColorManager";
import useGateLayers from "../hooks/useGateLayers";
import { useChartArrayGrid } from "../hooks/useChartArrayGrid";
import {
    useDensityGridCells,
    useDensityGridContours,
    useDensityGridViewState,
} from "../hooks/useDensityGridCells";
import { useChartID, useConfig, useFilteredIndices } from "../hooks";
import { useDeckSelectionMouseRebind } from "../hooks/useDeckSelectionMouseRebind";
import { useOuterContainer } from "../screen_state";
import type { DualContourLegacyConfig } from "../contour_state";
import type { ScatterPlotConfig2D } from "../scatter_state";
import { useSpatialLayers } from "../spatial_context";
import ChartArrayLayout from "./ChartArrayLayout";
import {
    cloneDeckLayer,
    getDensityGridViewId,
    getDensityGridViewStates,
    getSerializableViewState,
    shouldDrawLayerInDeckDensityGrid,
} from "./densityGridUtils";

type DensityGridConfig = ScatterPlotConfig2D & DualContourLegacyConfig;

export default function DeckDensityGridComponent() {
    const chartId = useChartID();
    const config = useConfig<DensityGridConfig>();
    const rows = useFilteredIndices();
    const outerContainer = useOuterContainer();
    const deckRef = useRef<any>(null);
    const {
        scatterProps: { scatterplotLayer, greyScatterplotLayer, setScatterKeyboardActive },
        selectionLayer,
    } = useSpatialLayers();
    useDeckSelectionMouseRebind(outerContainer, selectionLayer, deckRef, {
        canvasKey: "grid",
    });
    const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();
    const { cells, densityFields, configuredFieldCount } = useDensityGridCells();

    const grid = useChartArrayGrid({
        chartId,
        cells,
        getViewId: (id, cell, index) => getDensityGridViewId(id, cell.key, index),
    });

    const contourLayers = useDensityGridContours(densityFields, grid.visibleCellIndices);
    const referenceCellWidth = grid.metrics.cellBounds[0]?.width ?? 0;
    const referenceCellHeight = grid.metrics.cellBounds[0]?.height ?? 0;
    const radiusPixels = useDensityGridViewState(referenceCellWidth, referenceCellHeight);

    const views = useMemo(
        () =>
            grid.visibleCellIndices
                .map((index) => {
                    const bounds = grid.getCellBounds(index);
                    const viewId = grid.viewIds[index];
                    if (!bounds || bounds.width <= 0 || bounds.height <= 0 || !viewId) return null;
                    return new OrthographicView({
                        id: viewId,
                        x: bounds.x,
                        y: bounds.y,
                        width: bounds.width,
                        height: bounds.height,
                        flipY: false,
                        controller: { dragPan: controllerOptions.dragPan },
                    });
                })
                .filter((view): view is OrthographicView => view !== null),
        [grid.visibleCellIndices, grid.viewIds, grid.metrics.cellBounds, controllerOptions.dragPan],
    );

    const viewState = useMemo(
        () => getDensityGridViewStates(grid.visibleViewIds, config.viewState),
        [grid.visibleViewIds, config.viewState],
    );

    const layers = useMemo(
        () =>
            grid.visibleCellIndices.flatMap((index) => {
                const contourLayer = contourLayers[index];
                const viewId = grid.viewIds[index];
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
                return [greyLayer, scatterLayer];
            }),
        [
            grid.visibleCellIndices,
            grid.viewIds,
            contourLayers,
            radiusPixels,
            scatterplotLayer,
            greyScatterplotLayer,
        ],
    );

    const gateLayers = useMemo(
        () => [gateDisplayLayer, gateLabelLayer].filter((layer) => layer !== null),
        [gateDisplayLayer, gateLabelLayer],
    );

    const allLayers = useMemo(
        () =>
            selectionLayer
                ? [...layers, ...gateLayers, selectionLayer]
                : [...layers, ...gateLayers],
        [layers, gateLayers, selectionLayer],
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
                    title={field.name || field.field}
                >
                    {field.name || field.field}
                </div>
            );
        },
        [densityFields],
    );

    const hasCanvas =
        grid.hasCanvas && Number.isFinite(radiusPixels) && radiusPixels > 0;

    const deckOverlay = useMemo(
        () =>
            hasCanvas ? (
                <DeckGL
                    ref={deckRef}
                    controller={false}
                    layerFilter={({ layer, viewport }) => {
                        const viewId =
                            layer.props && "viewId" in layer.props && typeof layer.props.viewId === "string"
                                ? layer.props.viewId
                                : undefined;
                        return shouldDrawLayerInDeckDensityGrid(
                            { id: layer.id, props: viewId === undefined ? undefined : { viewId } },
                            viewport.id,
                        );
                    }}
                    layers={allLayers}
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
        [hasCanvas, allLayers, views, viewState, config],
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
                ref={grid.scrollContainerRef}
                className="min-h-0 min-w-0 flex-1 overflow-auto overflow-x-hidden"
                onMouseDown={() => setScatterKeyboardActive(true)}
                onMouseEnter={() => setScatterKeyboardActive(true)}
                onMouseLeave={() => setScatterKeyboardActive(false)}
            >
                <ChartArrayLayout
                    cellCount={grid.cellCount}
                    cellKeys={grid.cellKeys}
                    layoutRef={grid.layoutRef}
                    canvasOverlay={deckOverlay}
                    renderCell={renderCell}
                />
            </div>
        </div>
    );
}
