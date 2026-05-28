import DeckGL from "@deck.gl/react";
import { OrthographicView, type OrthographicViewState } from "@deck.gl/core";
import type { Device } from "@luma.gl/core";
import { action } from "mobx";
import { useCallback, useEffect, useLayoutEffect, useMemo, useRef } from "react";
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
import ChartArrayCellLabel from "./ChartArrayCellLabel";
import { ChartArrayEmptyState } from "./ChartArrayEmptyState";
import { buildChartArrayDeckLayers } from "./chartArrayDeckLayers";
import { tagDeckLayerViewportScope } from "./deckLayerViewportScope";
import {
    ChartArrayRenderCache,
    getChartArrayStaticCacheKey,
} from "@/webgl/chartArrayRenderCache";
import { heatmapContourBackend } from "@/webgl/contourBackend";
import {
    getDensityGridViewId,
    getDensityGridViewStates,
    shouldDrawLayerInDeckDensityGrid,
} from "./densityGridUtils";

type DensityGridConfig = ScatterPlotConfig2D & DualContourLegacyConfig;
const ENABLE_STATIC_GRID_BUFFER =
    (import.meta.env.VITE_MDV_DENSITY_GRID_STATIC_BUFFER ?? "true") !== "false";
const DENSITY_GRID_EMPTY_STATE_MESSAGES = {
    noCellsConfigured: "Choose density fields to build the grid.",
    loadingCells: "Loading density fields...",
} as const;

function applyDeckViewStateChange(
    update: OrthographicViewState | Record<string, OrthographicViewState>,
    current?: OrthographicViewState,
): OrthographicViewState {
    if (
        update &&
        typeof update === "object" &&
        "target" in update &&
        Array.isArray((update as OrthographicViewState).target)
    ) {
        const state = update as OrthographicViewState;
        return {
            target: state.target,
            zoom: state.zoom,
            minZoom: state.minZoom,
            maxZoom: state.maxZoom,
        };
    }
    const record = update as Record<string, OrthographicViewState>;
    const next = Object.values(record).find(
        (viewState) =>
            viewState &&
            Array.isArray(viewState.target) &&
            Number.isFinite(Number(viewState.zoom)),
    );
    if (next) {
        return {
            target: next.target,
            zoom: next.zoom,
            minZoom: next.minZoom,
            maxZoom: next.maxZoom,
        };
    }
    if (current) {
        return {
            target: current.target,
            zoom: current.zoom,
            minZoom: current.minZoom,
            maxZoom: current.maxZoom,
        };
    }
    return { target: [0, 0, 0], zoom: 0 };
}

function resolveDeckDevice(deck: unknown): Device | null {
    if (!deck || typeof deck !== "object") return null;
    const typed = deck as Record<string, unknown>;
    if (typed.device) return typed.device as Device;
    const animationLoop = typed.animationLoop as Record<string, unknown> | undefined;
    if (animationLoop?.device) return animationLoop.device as Device;
    const layerManager = typed.layerManager as Record<string, unknown> | undefined;
    const context = layerManager?.context as Record<string, unknown> | undefined;
    if (context?.device) return context.device as Device;
    return null;
}

export default function DeckDensityGridComponent() {
    const chartId = useChartID();
    const config = useConfig<DensityGridConfig>();
    const outerContainer = useOuterContainer();
    const deckRef = useRef<any>(null);
    const renderCacheRef = useRef<ChartArrayRenderCache | null>(null);
    const {
        scatterProps: { scatterplotLayer, greyScatterplotLayer, setScatterKeyboardActive },
        selectionLayer,
    } = useSpatialLayers();
    useDeckSelectionMouseRebind(outerContainer, selectionLayer, deckRef, {
        canvasKey: "grid",
    });
    const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();
    const { cells, densityFields, configuredFieldCount } = useDensityGridCells();

    useLayoutEffect(() => {
        if (configuredFieldCount === 0 && config.density_mode === "grid") {
            action(() => {
                config.density_mode = "overlay";
            })();
        }
    }, [configuredFieldCount, config]);

    const grid = useChartArrayGrid({
        chartId,
        cells,
        getViewId: (id, cell, index) => getDensityGridViewId(id, cell.key, index),
    });

    const contourLayers = useDensityGridContours(densityFields, grid.visibleCellIndices);
    const referenceCellWidth = grid.metrics.cellBounds[0]?.width ?? 0;
    const referenceCellHeight = grid.metrics.cellBounds[0]?.height ?? 0;
    const radiusPixels = useDensityGridViewState(referenceCellWidth, referenceCellHeight);

    const views = useMemo(() => {
        const built = grid.visibleCellIndices
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
            .filter((view): view is OrthographicView => view !== null);
        if (built.length > 0 || grid.cellCount === 0) return built;
        const viewId = grid.viewIds[0];
        const { width, height } = grid.metrics.rootSize;
        if (!viewId || width <= 0 || height <= 0) return built;
        return [
            new OrthographicView({
                id: viewId,
                x: 0,
                y: 0,
                width,
                height,
                flipY: false,
                controller: { dragPan: controllerOptions.dragPan },
            }),
        ];
    }, [
        grid.visibleCellIndices,
        grid.viewIds,
        grid.getCellBounds,
        grid.metrics.rootSize,
        grid.cellCount,
        controllerOptions.dragPan,
    ]);

    const viewState = useMemo(
        () => getDensityGridViewStates(grid.visibleViewIds, config.viewState),
        [grid.visibleViewIds, config.viewState],
    );

    const devicePixelRatio =
        typeof window !== "undefined" ? window.devicePixelRatio : 1;
    const staticCacheKey = useMemo(
        () =>
            getChartArrayStaticCacheKey({
                width: grid.metrics.rootSize.width,
                height: grid.metrics.rootSize.height,
                devicePixelRatio,
            }),
        [
            grid.metrics.rootSize.width,
            grid.metrics.rootSize.height,
            devicePixelRatio,
        ],
    );

    useEffect(() => {
        if (!ENABLE_STATIC_GRID_BUFFER) {
            return;
        }
        const deck = deckRef.current?.deck;
        const device = resolveDeckDevice(deck);
        if (!device) {
            return;
        }
        if (!renderCacheRef.current) {
            renderCacheRef.current = new ChartArrayRenderCache();
        }
        const cssWidth = Math.max(1, Math.round(grid.metrics.rootSize.width));
        const cssHeight = Math.max(1, Math.round(grid.metrics.rootSize.height));
        renderCacheRef.current.getOrCreateStaticPass({
            device,
            cacheKey: staticCacheKey,
            width: Math.max(1, Math.round(cssWidth * devicePixelRatio)),
            height: Math.max(1, Math.round(cssHeight * devicePixelRatio)),
            cssWidth,
            cssHeight,
        });
    }, [
        staticCacheKey,
        grid.metrics.rootSize.width,
        grid.metrics.rootSize.height,
        devicePixelRatio,
    ]);

    useEffect(
        () => () => {
            renderCacheRef.current?.dispose();
            renderCacheRef.current = null;
        },
        [],
    );

    const buildPerViewportLayer = useCallback(
        (index: number, viewId: string) => {
            const contourLayer = contourLayers[index];
            if (!contourLayer) return null;
            return tagDeckLayerViewportScope(
                heatmapContourBackend.buildDynamicLayer({
                    contourLayer,
                    id: `${viewId}-density`,
                    radiusPixels,
                    debounce: 250,
                    weightsTextureSize: 256,
                }),
                "per-viewport",
                { viewId },
            );
        },
        [contourLayers, radiusPixels],
    );

    const allLayers = useMemo(
        () =>
            buildChartArrayDeckLayers({
                visibleCellIndices: grid.visibleCellIndices,
                viewIds: grid.viewIds,
                greyScatterplotLayer: greyScatterplotLayer ?? null,
                scatterplotLayer: scatterplotLayer ?? null,
                chartSharedLayers: [gateDisplayLayer, gateLabelLayer, selectionLayer],
                buildPerViewportLayer,
            }),
        [
            grid.visibleCellIndices,
            grid.viewIds,
            greyScatterplotLayer,
            scatterplotLayer,
            gateDisplayLayer,
            gateLabelLayer,
            selectionLayer,
            buildPerViewportLayer,
        ],
    );

    const renderCell = useCallback(
        (index: number) => {
            const field = densityFields[index];
            if (!field) return null;
            const color = getFieldColor(field.field);
            return (
                <ChartArrayCellLabel
                    label={field.name || field.field}
                    accentRgb={color}
                />
            );
        },
        [densityFields],
    );

    const hasCanvas = grid.hasCanvas && views.length > 0;

    const deckOverlay = useMemo(
        () =>
            hasCanvas ? (
                <DeckGL
                    ref={deckRef}
                    controller={false}
                    layerFilter={({ layer, viewport }) =>
                        shouldDrawLayerInDeckDensityGrid(
                            { id: layer.id, props: layer.props as Record<string, unknown> },
                            viewport.id,
                        )
                    }
                    layers={allLayers}
                    views={views}
                    viewState={viewState}
                    useDevicePixels={true}
                    onViewStateChange={({
                        viewState: nextViewState,
                    }: {
                        viewState:
                            | OrthographicViewState
                            | Record<string, OrthographicViewState>;
                    }) => {
                        action(() => {
                            config.viewState = applyDeckViewStateChange(
                                nextViewState,
                                config.viewState,
                            );
                        })();
                    }}
                    getCursor={({ isDragging }) => (isDragging ? "grabbing" : "crosshair")}
                />
            ) : null,
        [hasCanvas, allLayers, views, viewState, config],
    );

    if (configuredFieldCount === 0) {
        return (
            <ChartArrayEmptyState
                configuredCellCount={0}
                loadedCellCount={0}
                noCellsConfiguredMessage={DENSITY_GRID_EMPTY_STATE_MESSAGES.noCellsConfigured}
                loadingCellsMessage={DENSITY_GRID_EMPTY_STATE_MESSAGES.loadingCells}
            />
        );
    }
    if (densityFields.length === 0) {
        return (
            <ChartArrayEmptyState
                configuredCellCount={configuredFieldCount}
                loadedCellCount={0}
                noCellsConfiguredMessage={DENSITY_GRID_EMPTY_STATE_MESSAGES.noCellsConfigured}
                loadingCellsMessage={DENSITY_GRID_EMPTY_STATE_MESSAGES.loadingCells}
            />
        );
    }
    return (
        <div className="relative flex h-full w-full min-w-0 flex-col">
            <div
                ref={grid.scrollContainerRef}
                className="relative min-h-0 min-w-0 flex-1 overflow-auto overflow-x-hidden"
                onMouseDown={() => setScatterKeyboardActive(true)}
                onMouseEnter={() => setScatterKeyboardActive(true)}
                onMouseLeave={() => setScatterKeyboardActive(false)}
            >
                <ChartArrayLayout
                    className="min-h-full w-full"
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
