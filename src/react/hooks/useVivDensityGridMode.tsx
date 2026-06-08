import { ColorPaletteExtension, DetailView } from "@hms-dbmi/viv";
import { action } from "mobx";
import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState, type ReactNode, type RefObject } from "react";
import { shallow } from "zustand/shallow";
import type { DualContourLegacyConfig } from "../contour_state";
import { getFieldColor } from "../fieldColorManager";
import useGateLayers from "./useGateLayers";
import { useChartArrayGrid } from "./useChartArrayGrid";
import {
    useDensityGridCells,
    useDensityGridContours,
} from "./useDensityGridCells";
import { useVivDensityGridViewState } from "./useVivDensityGridViewState";
import { useChartID, useConfig } from "../hooks";
import { useSpatialLayers } from "../spatial_context";
import ChartArrayLayout from "../components/ChartArrayLayout";
import ChartArrayCellLabel from "../components/ChartArrayCellLabel";
import { ChartArrayEmptyState } from "../components/ChartArrayEmptyState";
import {
    buildDeckLayerCacheKey,
    getDeckLayerTemplateSignature,
    useBuiltChartArrayDeckLayers,
} from "../components/chartArrayDeckLayerClones";
import { resolveDeckDevice } from "../components/deckDeviceUtils";
import { useDeckMountEpoch } from "./useDeckMountEpoch";
import { createHeatmapContourLayer } from "@/webgl/SpatialLayer";
import {
    tagDeckLayerPickOnlyInStaticComposite,
    tagDeckLayerViewportScope,
} from "../components/deckLayerViewportScope";
import {
    cloneDeckLayerForRender,
    DENSITY_GRID_EMPTY_STATE_MESSAGES,
    getVivGridDetailViewId,
} from "../components/densityGridUtils";
import { useJsonLayer } from "../components/vivJsonLayer";
import type { VivDensityGridViewerProps } from "../components/vivDensityGridViewerProps";
import type { VivRoiConfig } from "../components/VivMDVReact";
import { useLoader, type OME_TIFF, useChannelsStore, useViewerStore, useViewerStoreApi } from "../components/avivatorish/state";
import {
    getVivId,
    type ViewState as VivViewState,
    type ViewStates,
} from "../components/avivatorish/MDVivViewer";
import type { EditableGeoJsonLayer } from "@deck.gl-community/editable-layers";
import VivContrastExtension from "@/webgl/VivContrastExtension";
import { useOuterContainer } from "../screen_state";
import type { DeckGLProps, Deck, OrbitViewState, OrthographicViewState, PickingInfo } from "deck.gl";
import type { Layer } from "@deck.gl/core";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import { useOuterContainerDeckTooltip } from "./useOuterContainerDeckTooltip";
import {
    ChartArrayRenderCache,
    getChartArrayStaticCacheKey,
    getChartArrayStaticContentKey,
    type ChartArrayStaticPass,
} from "@/webgl/chartArrayRenderCache";
import { ChartArrayStaticPassRenderer } from "@/webgl/chartArrayStaticPassRenderer";
import {
    createChartArrayStaticCompositeLayer,
    getFramebufferCompositeBounds,
} from "@/webgl/chartArrayStaticBitmapLayers";
import { buildVivImageLayersForStaticPass } from "@/webgl/chartArrayStaticLayers";

function buildGridDeckLayersWithStaticComposite(
    allDeckLayers: Layer[],
    staticPass: ChartArrayStaticPass,
    staticContentFrame: number,
    viewState: VivViewState,
    viewStates: VivViewState[],
    visibleCellIndices: readonly number[],
    viewIds: readonly string[],
    rootSize: { width: number; height: number },
): Layer[] {
    const pickOnlyGeometryLayers = allDeckLayers
        .filter(isChartArrayScatterStaticLayer)
        .map((layer) => tagDeckLayerPickOnlyInStaticComposite(layer as CloneableDeckLayer));
    const baseLayers = allDeckLayers.filter((layer) => !isChartArrayScatterStaticLayer(layer));
    if (visibleCellIndices.length === 0) {
        return [...pickOnlyGeometryLayers, ...baseLayers];
    }
    const viewStateById = new Map(viewStates.map((vs) => [vs.id, vs]));
    const compositeLayers = visibleCellIndices.flatMap((index) => {
        const viewId = viewIds[index];
        if (!viewId) return [];
        const perViewState = viewStateById.get(viewId);
        const cellWidth = (perViewState as { width?: number } | undefined)?.width;
        const cellHeight = (perViewState as { height?: number } | undefined)?.height;
        const boundsWidth = typeof cellWidth === "number" && cellWidth > 0 ? cellWidth : rootSize.width;
        const boundsHeight =
            typeof cellHeight === "number" && cellHeight > 0 ? cellHeight : rootSize.height;
        if (boundsWidth <= 0 || boundsHeight <= 0) return [];
        const bounds = getFramebufferCompositeBounds(
            (perViewState ?? viewState) as OrthographicViewState,
            boundsWidth,
            boundsHeight,
        );
        return [
            createChartArrayStaticCompositeLayer(
                `chart-array-static-composite-${viewId}-${staticContentFrame}`,
                viewId,
                staticPass.colorTexture,
                bounds,
            ),
        ];
    });
    return [...compositeLayers, ...pickOnlyGeometryLayers, ...baseLayers];
}

export const VIV_SCATTER_DECK_KEY = "viv-scatter-deck";

const ENABLE_STATIC_GRID_BUFFER =
    (import.meta.env.VITE_MDV_DENSITY_GRID_STATIC_BUFFER ?? "true") !== "false";

function isChartArrayScatterStaticLayer(layer: Layer): boolean {
    return layer.id.startsWith("scatter_") || layer.id.startsWith("scatter-grey_");
}

/** Shared pan/zoom for the offscreen pass (viewport size is on the OrthographicView). */
function buildStaticPassViewState(viewState: VivViewState): OrthographicViewState {
    return {
        target: viewState.target,
        zoom: viewState.zoom,
        minZoom: viewState.minZoom,
        maxZoom: viewState.maxZoom,
    };
}

type CloneableDeckLayer = Layer & {
    clone: (props: Record<string, unknown>) => Layer;
};

function cloneLayersForStaticPass(layers: Layer[]): Layer[] {
    return layers.map((layer) => {
        const cloneable = layer as CloneableDeckLayer;
        if (!cloneable.clone) return layer;
        return cloneDeckLayerForRender(cloneable, {
            id: `${layer.id}-static-fbo`,
        });
    });
}

export type { VivDensityGridViewerProps } from "../components/vivDensityGridViewerProps";

export type VivDensityGridModeResult = {
    blocking: ReactNode | null;
    layout: ReactNode | null;
    viewer: VivDensityGridViewerProps | null;
    tooltipPortal: ReactNode;
    containerHandlers: {
        onPointerDown: () => void;
        onMouseDown: () => void;
        onMouseEnter: () => void;
        onMouseLeave: () => void;
    };
};

export function useVivDensityGridMode(
    enabled: boolean,
    deckContainerRef: RefObject<HTMLDivElement | null>,
): VivDensityGridModeResult | null {
    const chartId = useChartID();
    const ome = useLoader() as OME_TIFF["data"];
    const viewerStore = useViewerStoreApi();
    const viewState = useViewerStore((store) => store.viewState);
    const outerContainer = useOuterContainer();

    const {
        scatterProps: {
            scatterplotLayer,
            greyScatterplotLayer,
            getTooltip,
            setScatterKeyboardActive,
        },
        selectionLayer,
    } = useSpatialLayers();
    const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();
    const config = useConfig<VivRoiConfig & DualContourLegacyConfig>();
    const { showJson } = config;
    const jsonLayer = useJsonLayer(showJson, "chart-array");
    const { cells, densityFields, configuredFieldCount } = useDensityGridCells();

    useLayoutEffect(() => {
        if (!enabled) return;
        if (configuredFieldCount === 0 && config.density_mode === "grid") {
            action(() => {
                config.density_mode = "overlay";
            })();
        }
    }, [enabled, configuredFieldCount, config]);

    const deckMountEpoch = useDeckMountEpoch(enabled);

    const grid = useChartArrayGrid({
        chartId,
        cells: enabled ? cells : [],
        getViewId: (id, cell, index) => getVivGridDetailViewId(id, cell.key, index),
    });

    const layoutCellIndices = useMemo(
        () => (enabled ? Array.from({ length: grid.cellCount }, (_, index) => index) : []),
        [enabled, grid.cellCount],
    );

    const contourLayers = useDensityGridContours(
        enabled ? densityFields : [],
        grid.visibleCellIndices,
    );
    const referenceCellWidth = grid.metrics.cellBounds[0]?.width ?? 0;
    const referenceCellHeight = grid.metrics.cellBounds[0]?.height ?? 0;

    const radiusPixels = useVivDensityGridViewState(referenceCellWidth, referenceCellHeight);

    const { colors, contrastLimits, channelsVisible, selections, brightness, contrast } = useChannelsStore(
        ({ colors, contrastLimits, channelsVisible, selections, brightness, contrast }) => ({
            colors,
            contrastLimits,
            channelsVisible,
            selections,
            brightness,
            contrast,
        }),
        shallow,
    );

    const extensions = useMemo(() => [new ColorPaletteExtension(), new VivContrastExtension()], []);
    const layerConfig = useMemo(
        () => ({
            loader: ome,
            selections,
            contrastLimits,
            extensions,
            colors,
            channelsVisible,
            brightness,
            contrast,
        }),
        [ome, selections, contrastLimits, extensions, colors, channelsVisible, brightness, contrast],
    );

    const detailViews = useMemo(
        () =>
            layoutCellIndices
                .map((index) => {
                    const bounds = grid.getCellBounds(index);
                    const detailId = grid.viewIds[index];
                    if (!bounds || bounds.width <= 0 || bounds.height <= 0 || !detailId) return null;
                    return new DetailView({
                        id: detailId,
                        x: bounds.x,
                        y: bounds.y,
                        width: bounds.width,
                        height: bounds.height,
                        // @ts-expect-error viv runtime supports this, types do not
                        snapScaleBar: false,
                    });
                })
                .filter((view): view is DetailView => view !== null),
        [layoutCellIndices, grid.viewIds, grid.getCellBounds],
    );

    const viewStates = useMemo((): ViewStates => {
        if (!viewState) return [];
        return layoutCellIndices
            .map((index) => {
                const detailId = grid.viewIds[index];
                const bounds = grid.getCellBounds(index);
                if (!detailId) return null;
                return {
                    ...viewState,
                    id: detailId,
                    ...(bounds && bounds.width > 0 && bounds.height > 0
                        ? { width: bounds.width, height: bounds.height }
                        : {}),
                } satisfies VivViewState;
            })
            .filter((vs): vs is NonNullable<typeof vs> => vs !== null) satisfies ViewStates;
    }, [layoutCellIndices, grid.viewIds, grid.getCellBounds, viewState]);

    const buildPerViewportLayer = useCallback(
        (index: number, viewId: string) => {
            const contourLayer = contourLayers[index];
            if (!contourLayer) return null;
            return tagDeckLayerViewportScope(
                createHeatmapContourLayer({
                    ...contourLayer,
                    radiusPixels,
                    debounce: 250,
                    weightsTextureSize: 256,
                    id: `density_${getVivId(viewId)}`,
                }),
                "per-viewport",
                { viewId },
            );
        },
        [contourLayers, radiusPixels],
    );

    const chartArrayLayerBuildOptions = useMemo(
        () => ({
            visibleCellIndices: grid.visibleCellIndices,
            viewIds: grid.viewIds,
            greyScatterplotLayer: greyScatterplotLayer ?? null,
            scatterplotLayer: scatterplotLayer ?? null,
            chartSharedLayers: [jsonLayer, gateDisplayLayer, gateLabelLayer, selectionLayer],
            buildPerViewportLayer,
        }),
        [
            grid.visibleCellIndices,
            grid.viewIds,
            greyScatterplotLayer,
            scatterplotLayer,
            jsonLayer,
            gateDisplayLayer,
            gateLabelLayer,
            selectionLayer,
            buildPerViewportLayer,
        ],
    );

    const chartArrayDeckLayersCacheKey = useMemo(
        () =>
            buildDeckLayerCacheKey([
                `epoch=${deckMountEpoch}`,
                grid.visibleCellIndices.join(","),
                grid.viewIds.join(","),
                String(radiusPixels),
                greyScatterplotLayer
                    ? getDeckLayerTemplateSignature(greyScatterplotLayer)
                    : "",
                scatterplotLayer ? getDeckLayerTemplateSignature(scatterplotLayer) : "",
                jsonLayer ? getDeckLayerTemplateSignature(jsonLayer) : "",
                gateDisplayLayer ? getDeckLayerTemplateSignature(gateDisplayLayer) : "",
                gateLabelLayer ? getDeckLayerTemplateSignature(gateLabelLayer) : "",
                selectionLayer ? getDeckLayerTemplateSignature(selectionLayer) : "",
                ...grid.visibleCellIndices.map(
                    (index) => densityFields[index]?.field ?? `cell-${index}`,
                ),
            ]),
        [
            deckMountEpoch,
            grid.visibleCellIndices,
            grid.viewIds,
            radiusPixels,
            greyScatterplotLayer,
            scatterplotLayer,
            jsonLayer,
            gateDisplayLayer,
            gateLabelLayer,
            selectionLayer,
            densityFields,
        ],
    );

    const allDeckLayers = useBuiltChartArrayDeckLayers(
        chartArrayLayerBuildOptions,
        chartArrayDeckLayersCacheKey,
        deckMountEpoch,
    );

    const gridSelectionLayer = useMemo(
        () =>
            allDeckLayers.find((layer) => layer.id.startsWith("selection_")) as
                | EditableGeoJsonLayer
                | undefined,
        [allDeckLayers],
    );

    const renderCacheRef = useRef<ChartArrayRenderCache | null>(null);
    const staticPassRendererRef = useRef<ChartArrayStaticPassRenderer | null>(null);
    const deckInstanceRef = useRef<Deck | null>(null);
    const [staticPassRevision, setStaticPassRevision] = useState(0);
    const [staticPassFailed, setStaticPassFailed] = useState(false);
    const staticPassRef = useRef<{
        pass: ChartArrayStaticPass;
        renderKey: string;
        contentFrame: number;
    } | null>(null);
    const staticContentFrameRef = useRef(0);
    const devicePixelRatio =
        typeof window !== "undefined" ? window.devicePixelRatio : 1;

    // Match reference cell 0 used by useVivDensityGridViewState for zoom fit.
    const staticPassPixelSize = useMemo(() => {
        if (grid.cellCount === 0 || referenceCellWidth <= 0 || referenceCellHeight <= 0) {
            return { width: 0, height: 0 };
        }
        return {
            width: Math.max(1, Math.round(referenceCellWidth)),
            height: Math.max(1, Math.round(referenceCellHeight)),
        };
    }, [grid.cellCount, referenceCellWidth, referenceCellHeight]);

    const vivLayersForStaticPass = useMemo(
        () =>
            ome && staticPassPixelSize.width > 0 && staticPassPixelSize.height > 0
                ? buildVivImageLayersForStaticPass({
                      width: staticPassPixelSize.width,
                      height: staticPassPixelSize.height,
                      layerConfig,
                  })
                : [],
        [ome, staticPassPixelSize.width, staticPassPixelSize.height, layerConfig],
    );

    const staticSourceLayers = useMemo(
        () => [
            ...vivLayersForStaticPass,
            ...allDeckLayers.filter(isChartArrayScatterStaticLayer),
        ],
        [allDeckLayers, vivLayersForStaticPass],
    );

    const vivStaticLayerKey = useMemo(
        () =>
            buildDeckLayerCacheKey([
                String(staticPassPixelSize.width),
                String(staticPassPixelSize.height),
                JSON.stringify(layerConfig.channelsVisible),
                JSON.stringify(layerConfig.contrastLimits),
                JSON.stringify(layerConfig.selections),
                JSON.stringify(layerConfig.colors),
                String(layerConfig.brightness),
                String(layerConfig.contrast),
            ]),
        [staticPassPixelSize.width, staticPassPixelSize.height, layerConfig],
    );

    const staticLayerIds = useMemo(
        () => staticSourceLayers.map((layer) => layer.id),
        [staticSourceLayers],
    );

    const staticBufferKey = useMemo(
        () =>
            getChartArrayStaticCacheKey({
                width: staticPassPixelSize.width,
                height: staticPassPixelSize.height,
                devicePixelRatio,
            }),
        [staticPassPixelSize.width, staticPassPixelSize.height, devicePixelRatio],
    );

    const staticContentKey = useMemo(
        () => (viewState ? getChartArrayStaticContentKey(viewState, staticLayerIds) : null),
        [viewState, staticLayerIds],
    );

    const staticRenderKey = useMemo(() => {
        if (!staticContentKey) return null;
        const viewIdsKey = grid.viewIds.join(",");
        const rootSizeKey = `${grid.metrics.rootSize.width}x${grid.metrics.rootSize.height}`;
        const visibleCellIndicesKey = grid.visibleCellIndices.join(",");
        return `${staticBufferKey}|${staticContentKey}|${chartArrayDeckLayersCacheKey}|${vivStaticLayerKey}|${visibleCellIndicesKey}|${viewIdsKey}|${rootSizeKey}`;
    }, [
        staticBufferKey,
        staticContentKey,
        chartArrayDeckLayersCacheKey,
        vivStaticLayerKey,
        grid.visibleCellIndices,
        grid.viewIds,
        grid.metrics.rootSize,
    ]);

    const onDeckInstance = useCallback((deck: Deck | null) => {
        deckInstanceRef.current = deck;
        if (!deck) {
            renderCacheRef.current?.dispose();
            renderCacheRef.current = null;
            staticPassRendererRef.current?.finalize();
            staticPassRendererRef.current = null;
            staticPassRef.current = null;
            return;
        }
    }, []);

    // Rasterize shared geometry once per staticRenderKey; ref + revision avoids setState loops.
    useLayoutEffect(() => {
        if (!enabled || !ENABLE_STATIC_GRID_BUFFER || staticPassFailed || !staticRenderKey) {
            if (staticPassRef.current !== null) {
                staticPassRef.current = null;
                setStaticPassRevision((revision) => revision + 1);
            }
            return;
        }
        if (staticPassRef.current?.renderKey === staticRenderKey) {
            return;
        }
        const deck = deckInstanceRef.current;
        if (!deck || !viewState || staticSourceLayers.length === 0) {
            return;
        }
        const device = resolveDeckDevice(deck);
        if (!device) {
            return;
        }

        const cssWidth = staticPassPixelSize.width;
        const cssHeight = staticPassPixelSize.height;
        if (cssWidth <= 0 || cssHeight <= 0) {
            return;
        }
        const physicalWidth = Math.max(1, Math.round(cssWidth * devicePixelRatio));
        const physicalHeight = Math.max(1, Math.round(cssHeight * devicePixelRatio));
        const staticViewState = buildStaticPassViewState(viewState);
        const layersForStaticPass = cloneLayersForStaticPass(staticSourceLayers);

        if (!renderCacheRef.current) {
            renderCacheRef.current = new ChartArrayRenderCache();
        }
        const pass = renderCacheRef.current.getOrCreateStaticPass({
            device,
            cacheKey: staticBufferKey,
            width: physicalWidth,
            height: physicalHeight,
            cssWidth,
            cssHeight,
        });

        try {
            if (!staticPassRendererRef.current) {
                staticPassRendererRef.current = new ChartArrayStaticPassRenderer();
            }
            staticPassRendererRef.current.render({
                device,
                target: pass.framebuffer,
                cssWidth,
                cssHeight,
                viewState: staticViewState,
                layers: layersForStaticPass,
            });
            staticContentFrameRef.current += 1;
            staticPassRef.current = {
                pass,
                renderKey: staticRenderKey,
                contentFrame: staticContentFrameRef.current,
            };
            setStaticPassRevision((revision) => revision + 1);
        } catch {
            renderCacheRef.current?.dispose();
            renderCacheRef.current = null;
            staticPassRendererRef.current?.finalize();
            staticPassRendererRef.current = null;
            staticPassRef.current = null;
            setStaticPassFailed(true);
        }
    }, [
        enabled,
        staticPassFailed,
        staticRenderKey,
        staticBufferKey,
        staticSourceLayers,
        viewState,
        staticPassPixelSize.width,
        staticPassPixelSize.height,
        devicePixelRatio,
    ]);

    const staticCompositeActive = useMemo(() => {
        staticPassRevision;
        if (!ENABLE_STATIC_GRID_BUFFER || staticPassFailed || !enabled) {
            return false;
        }
        const cached = staticPassRef.current;
        return Boolean(cached && cached.renderKey === staticRenderKey);
    }, [staticPassRevision, staticPassFailed, staticRenderKey, enabled]);

    const deckLayersForRender = useMemo(() => {
        // Invalidation only: re-run when the offscreen raster pass finishes (staticPassRef is a ref).
        staticPassRevision;
        if (!ENABLE_STATIC_GRID_BUFFER || staticPassFailed || !enabled || !viewState) {
            return allDeckLayers;
        }
        const cached = staticPassRef.current;
        if (!cached || cached.renderKey !== staticRenderKey) {
            return allDeckLayers;
        }
        return buildGridDeckLayersWithStaticComposite(
            allDeckLayers,
            cached.pass,
            cached.contentFrame,
            viewState,
            viewStates,
            grid.visibleCellIndices,
            grid.viewIds,
            grid.metrics.rootSize,
        );
    }, [
        allDeckLayers,
        enabled,
        staticPassFailed,
        staticPassRevision,
        staticRenderKey,
        viewState,
        viewStates,
        grid.visibleCellIndices,
        grid.viewIds,
        grid.metrics.rootSize,
    ]);

    useEffect(
        () => () => {
            renderCacheRef.current?.dispose();
            staticPassRendererRef.current?.finalize();
            renderCacheRef.current = null;
            staticPassRendererRef.current = null;
        },
        [],
    );

    const getTooltipContent = useCallback(
        (info: PickingInfo) =>
            getCombinedScatterTooltip(info, {
                gateDisplayLayerId: gateDisplayLayer?.id,
                gateLabelLayerId: gateLabelLayer?.id,
                getPointTooltip: getTooltip,
            }),
        [gateDisplayLayer?.id, gateLabelLayer?.id, getTooltip],
    );

    const {
        clearTooltip,
        getTooltip: getPortalTooltip,
        suppressTooltipUntilPointerUp,
        tooltipPortal,
    } = useOuterContainerDeckTooltip(getTooltipContent, deckContainerRef);

    const deckProps: Partial<DeckGLProps> = useMemo(
        () => ({
            getTooltip: getPortalTooltip,
            layers: deckLayersForRender,
            controller: {
                doubleClickZoom: false,
                dragPan: controllerOptions.dragPan,
            },
        }),
        [deckLayersForRender, getPortalTooltip, controllerOptions.dragPan],
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

    const canRenderDeck =
        enabled &&
        grid.hasCanvas &&
        !!ome &&
        !!viewState &&
        Number.isFinite(radiusPixels) &&
        radiusPixels > 0 &&
        grid.cellCount > 0;

    const onViewStateChange = useCallback(
        (e: { viewId: string; viewState: OrthographicViewState | OrbitViewState }) => {
            viewerStore.setState({
                viewState: { ...e.viewState, id: e.viewId },
            });
        },
        [viewerStore],
    );

    const containerHandlers = useMemo(
        () => ({
            onPointerDown: suppressTooltipUntilPointerUp,
            onMouseDown: () => setScatterKeyboardActive(true),
            onMouseEnter: () => setScatterKeyboardActive(true),
            onMouseLeave: () => {
                clearTooltip();
                setScatterKeyboardActive(false);
            },
        }),
        [suppressTooltipUntilPointerUp, setScatterKeyboardActive, clearTooltip],
    );

    return useMemo((): VivDensityGridModeResult | null => {
        if (!enabled) return null;

        const emptyMessages = DENSITY_GRID_EMPTY_STATE_MESSAGES;

        if (configuredFieldCount === 0) {
            return {
                blocking: (
                    <ChartArrayEmptyState
                        configuredCellCount={0}
                        loadedCellCount={0}
                        noCellsConfiguredMessage={emptyMessages.noCellsConfigured}
                        loadingCellsMessage={emptyMessages.loadingCells}
                    />
                ),
                layout: null,
                viewer: null,
                tooltipPortal: null,
                containerHandlers,
            };
        }
        if (!viewState) {
            return {
                blocking: (
                    <ChartArrayEmptyState
                        waitingForViewState={true}
                        configuredCellCount={configuredFieldCount}
                        loadedCellCount={densityFields.length}
                        noCellsConfiguredMessage={emptyMessages.noCellsConfigured}
                        loadingCellsMessage={emptyMessages.loadingCells}
                    />
                ),
                layout: null,
                viewer: null,
                tooltipPortal: null,
                containerHandlers,
            };
        }
        if (densityFields.length === 0) {
            return {
                blocking: (
                    <ChartArrayEmptyState
                        configuredCellCount={configuredFieldCount}
                        loadedCellCount={0}
                        noCellsConfiguredMessage={emptyMessages.noCellsConfigured}
                        loadingCellsMessage={emptyMessages.loadingCells}
                    />
                ),
                layout: null,
                viewer: null,
                tooltipPortal: null,
                containerHandlers,
            };
        }

        const layout = (
            <div
                ref={grid.scrollContainerRef}
                className="relative min-h-0 min-w-0 flex-1 overflow-auto overflow-x-hidden"
            >
                <ChartArrayLayout
                    className="min-h-full w-full"
                    cellCount={grid.cellCount}
                    cellKeys={grid.cellKeys}
                    layoutRef={grid.layoutRef}
                    canvasOverlay={null}
                    renderCell={renderCell}
                />
            </div>
        );

        const viewer: VivDensityGridViewerProps | null = canRenderDeck
            ? {
                  outerContainer,
                  selectionLayer: gridSelectionLayer ?? selectionLayer,
                  views: detailViews,
                  layerProps: detailViews.map(() => layerConfig),
                  viewStates,
                  onViewStateChange,
                  deckProps,
                  onDeckInstance,
                  vivImageLayersPickOnlyInStaticComposite: staticCompositeActive,
              }
            : null;

        return {
            blocking: null,
            layout,
            viewer,
            tooltipPortal,
            containerHandlers,
        };
    }, [
        enabled,
        configuredFieldCount,
        densityFields.length,
        viewState,
        canRenderDeck,
        grid.scrollContainerRef,
        grid.cellCount,
        grid.cellKeys,
        grid.layoutRef,
        renderCell,
        outerContainer,
        gridSelectionLayer,
        selectionLayer,
        detailViews,
        layerConfig,
        viewStates,
        onViewStateChange,
        deckProps,
        onDeckInstance,
        staticCompositeActive,
        tooltipPortal,
        containerHandlers,
    ]);
}
