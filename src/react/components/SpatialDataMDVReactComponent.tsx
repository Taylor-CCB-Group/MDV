import { ColorPaletteExtension } from "@hms-dbmi/viv";
import { Button } from "@mui/material";
import { viewStateFromBounds } from "@spatialdata/core";
import { AlphaFromIntensityExtension, findRenderStackEntry, getRenderStackEntryIds } from "@spatialdata/layers";
import type { EntryStatus, RenderStackLayerInputs } from "@spatialdata/vis";
import {
    resolvePickedEntryId,
    type ViewState as SpatialCanvasViewState,
    type SpatialFeaturePickEvent,
    SpatialViewer,
    useSpatialCanvasRendererFromLayerInputs,
    type VivImageLayerContext,
    type VivImagePropsResolver,
} from "@spatialdata/vis";
import type { DeckGLProps, DeckGLRef, Layer, OrthographicViewState, PickingInfo } from "deck.gl";
import { observer } from "mobx-react-lite";
import { Profiler, useCallback, useEffect, useMemo, useRef, useState } from "react";
import { createPortal } from "react-dom";
import { ErrorBoundary, type FallbackProps } from "react-error-boundary";

import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import { GatedOrthographicController } from "@/lib/viewGestureGate";
import {
    type DeckErrorHandlers,
    resolveDeckBackend,
    useDeckDevice,
    useDeckDeviceProps,
} from "@/react/spatialdata/deck_backend";
import { ensureParquetWorker } from "@/react/spatialdata/ensureParquetWorker";
import { entryStatusForDisplay } from "@/react/spatialdata/entry_status";
import { useProjectMdvFieldSpecs } from "@/react/spatialdata/field_spec_projection";
import { deckHostLayerId } from "@/react/spatialdata/host_overlay_ids";
import { createImageLayerRegistry } from "@/react/spatialdata/image_layer_registry";
import { applyIsolatedGroupFix } from "@/react/spatialdata/isolated_group_fix";
import { LayerFailureTracker } from "@/react/spatialdata/layer_failures";
import { onSpatialProfilerRender } from "@/react/spatialdata/perf";
import { createPointsLayerRegistry } from "@/react/spatialdata/points_layer_registry";
import { describeGpu } from "@/react/spatialdata/render_report";
import {
    createMdvHostLayerResolver,
    type MdvDeckOverlayLayers,
    useRenderStackAdapter,
    useRenderStackStructure,
} from "@/react/spatialdata/render_stack_adapter";
import { seedRenderStackFromSpatialData } from "@/react/spatialdata/render_stack_control";
import { renderStackEntryDisplayName } from "@/react/spatialdata/render_stack_display";
import { useChartSpatialData } from "@/react/spatialdata/shared_spatial_data";
import type { SpatialDataState } from "@/react/spatialdata/spatial_data_reader";
import { formatSpatialFeatureTooltipHtml } from "@/react/spatialdata/spatial_feature_tooltip";
import { useAssociatedShapesLayerInputs } from "@/react/spatialdata/table_association";
import { toMdvViewState, toSpatialViewState } from "@/react/spatialdata/view_state_bridge";
import { imageAlphaFromIntensity, vivImageLayerProps } from "@/react/spatialdata/viv_image_props";
import VivContrastExtension from "@/webgl/VivContrastExtension";
import { useViewStateLink } from "../chartLinkHooks";
import { useChart } from "../context";
import type { DualContourLegacyConfig, HoveredContour } from "../contour_state";
import { useCategoryContourLegend, useFieldContourLegend } from "../contour_state";
import { useChartID, useChartSize, useConfig, useRegion } from "../hooks";
import { useDensityGridCells } from "../hooks/useDensityGridCells";
import useGateLayers from "../hooks/useGateLayers";
import { useOuterContainerDeckTooltip } from "../hooks/useOuterContainerDeckTooltip";
import { useRebindDeckEvents } from "../hooks/useRebindDeckEvents";
import { useOuterContainer } from "../screen_state";
import { SpatialAnnotationProvider, useSpatialLayers } from "../spatial_context";
import { useViewerStore, useViewerStoreApi, VivProvider } from "./avivatorish/state";
import { CrosshairArea } from "./CrosshairArea";
import { isDensityGridShown } from "./densityGridUtils";
import ErrorComponentReactWrapper from "./ErrorComponentReactWrapper";
import FieldContourLegend from "./FieldContourLegend";
import { isCrosshairShown } from "./pointerCrosshair";
import { isRegionGridShown } from "./regionGridUtils";
import SelectionOverlay from "./SelectionOverlay";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "./SpatialDataMDVReact";
import SpatialDensityGrid, { type SpatialDensityGridHost } from "./SpatialDensityGrid";
import SpatialRegionGrid, { type SpatialRegionGridHost } from "./SpatialRegionGrid";
import ViewGestureHint from "./ViewGestureHint";

// before any isolated group draws, in the single view or the density grid
applyIsolatedGroupFix();

type SpatialRegionMetadata = {
    spatial?: {
        file?: string;
        /** Where the store is served from, when not the project's `spatial/<file>`. */
        url?: string;
        coordinate_system?: string;
    };
};

function getSpatialRegionMetadata(region: unknown): SpatialRegionMetadata | null {
    if (!region || typeof region !== "object") return null;
    return region;
}

/**
 * Viv's image extensions and per-layer props for the chart's render stack, for the single view
 * and the grids alike.
 */
function useSpatialVivPassthrough(stack: SpatialDataMdvReactConfig["renderStack"]) {
    const extensions = useMemo(() => [new ColorPaletteExtension(), new VivContrastExtension()], []);
    // For image layers with "opacity from brightness" on, and only those (see
    // `imageAlphaFromIntensity`). Both lists keep their identity, so deck sees no change of
    // extensions until the option is turned on or off.
    const extensionsWithAlpha = useMemo(() => [...extensions, new AlphaFromIntensityExtension()], [extensions]);

    const vivImagePropsResolver = useCallback<VivImagePropsResolver>(
        (ctx: VivImageLayerContext) => vivImageLayerProps(stack, ctx.layerId),
        [stack],
    );

    return useMemo(
        () => ({
            vivImageExtensions: extensions,
            vivImageExtensionResolver: (ctx: VivImageLayerContext) =>
                imageAlphaFromIntensity(stack, ctx.layerId) ? extensionsWithAlpha : extensions,
            vivImagePropsResolver,
        }),
        [extensions, extensionsWithAlpha, stack, vivImagePropsResolver],
    );
}

const SpatialCanvasFromRenderStack = observer(function SpatialCanvasFromRenderStack({
    spatialData,
    coordinateSystem,
    spatialViewState,
    onSpatialViewStateChange,
    deckLayers: hostDeckLayers,
    layers,
    layerOrder,
    width,
    height,
    deckProps,
    onFeatureHover,
    onFeatureHoverEnd,
    densityGrid,
}: {
    spatialData: NonNullable<SpatialDataState["spatialData"]>;
    coordinateSystem: string;
    spatialViewState: ReturnType<typeof toSpatialViewState>;
    onSpatialViewStateChange: (next: SpatialCanvasViewState) => void;
    deckLayers: ReturnType<typeof useRenderStackAdapter>["deckLayers"];
    layers: ReturnType<typeof useRenderStackAdapter>["layers"];
    layerOrder: ReturnType<typeof useRenderStackAdapter>["layerOrder"];
    width: number;
    height: number;
    deckProps: Partial<DeckGLProps>;
    /** The feature's pick, and what its `x` / `y` are relative to when not the chart's canvas. */
    onFeatureHover: (event: SpatialFeaturePickEvent, anchor?: HTMLElement) => void;
    onFeatureHoverEnd: () => void;
    /** Set to draw each density field in its own cell instead of one canvas. */
    densityGrid?: SpatialDensityGridHost;
}) {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const stack = config.renderStack;
    void chart.renderStackGeneration;
    const stackStructure = useRenderStackStructure(stack);

    const deckRef = useRef<DeckGLRef | null>(null);
    useRebindDeckEvents(deckRef);

    const vivPassthrough = useSpatialVivPassthrough(stack);

    // What sd.js reports about each layer, for the chart's status and the layer dialog, which
    // live outside this tree. Called from an effect, with each layer's status as it changes.
    const onEntryStatusChange = useCallback(
        (layerId: string, status: EntryStatus | undefined) => {
            const entry = stack && findRenderStackEntry(stack, layerId);
            const elementType = entry?.kind === "spatial" ? entry.source.elementType : undefined;
            chart.setEntryStatus(layerId, entryStatusForDisplay(status, coordinateSystem, elementType));
        },
        [chart, stack, coordinateSystem],
    );
    useEffect(() => () => chart.clearEntryStatuses(), [chart]);

    const renderer = useSpatialCanvasRendererFromLayerInputs({
        spatialData,
        coordinateSystem,
        onEntryStatusChange,
        layerInputs: { layers, layerOrder },
        renderOrder: layerOrder,
        viewState: spatialViewState,
        onViewStateChange: onSpatialViewStateChange,
        width,
        height,
        hostDeckLayers,
        sortDeckLayers: true,
        // Auto-fit is driven by MDV below, not the library's built-in one — see the
        // effect for why (the built-in can commit a degenerate view before bounds load).
        autoFit: false,
        vivPassthrough,
    });

    // Fit the view to the image extent when no viewState is set. We drive this
    // ourselves instead of the library's `autoFit`: the built-in fires the moment
    // `isBlocking` clears and, if the visible layers' world bounds aren't populated
    // at that instant, commits a degenerate `{target:[0,0], zoom:0}` that then sticks
    // (viewState is no longer null, so it never re-fits). Gating on `hasLayersDrawn`
    // guarantees the bounds are available before we commit; while they aren't we
    // leave viewState null so this retries on the re-renders that happen as data loads.
    const { getWorldBoundsForVisibleLayers, hasEnabledLayers, hasLayersDrawn, isBlocking } = renderer;
    useEffect(() => {
        if (spatialViewState !== null) return;
        if (!hasEnabledLayers || isBlocking || !hasLayersDrawn) return;
        if (width <= 0 || height <= 0) return;
        const bounds = getWorldBoundsForVisibleLayers();
        if (!bounds) return;
        onSpatialViewStateChange(viewStateFromBounds(bounds, width, height));
    }, [
        spatialViewState,
        hasEnabledLayers,
        hasLayersDrawn,
        isBlocking,
        getWorldBoundsForVisibleLayers,
        width,
        height,
        onSpatialViewStateChange,
    ]);

    useEffect(() => {
        chart.renderStackGeneration; //explicitly add dependency for side-effects
        if (!stack) {
            chart.setImageLayerRegistry(undefined);
            return;
        }
        chart.setImageLayerRegistry(
            createImageLayerRegistry(stack, renderer.getImageLoadedDataByElementKey, renderer.getLayerLoadState),
        );
        return () => chart.setImageLayerRegistry(undefined);
    }, [
        chart,
        stack,
        chart.renderStackGeneration,
        renderer.getImageLoadedDataByElementKey,
        renderer.getLayerLoadState,
    ]);

    // The points engine, published for the layer dialog — a separate React tree, so
    // it cannot reach the renderer result. No render-stack dependency: the engine is
    // stable for the life of the renderer and owns its own notifications, so this
    // publishes once rather than on every generation bump.
    useEffect(() => {
        chart.setRetryEntry(renderer.retryEntry);
        return () => chart.setRetryEntry(undefined);
    }, [chart, renderer.retryEntry]);

    useEffect(() => {
        chart.setPointsLayerRegistry(createPointsLayerRegistry(renderer.pointsEngine, renderer.resolvePointsTarget));
        return () => chart.setPointsLayerRegistry(undefined);
    }, [chart, renderer.pointsEngine, renderer.resolvePointsTarget]);

    const handleHover = useCallback(
        (info: PickingInfo, anchor?: HTMLElement) => {
            if (!info.picked || typeof info.x !== "number" || typeof info.y !== "number") {
                // Returning without telling the host means an unpicked hover leaves the
                // last tooltip on screen: deck reports moving OFF a feature this way, and
                // the only other clear is `onMouseLeave` on the container, so the stale
                // tooltip would follow the pointer around the canvas until it left.
                onFeatureHoverEnd();
                return;
            }
            // Looks through isolated groups to the member, and drops the view token.
            const entryId = resolvePickedEntryId(info);
            const event = renderer.getFeaturePickEvent(entryId, {
                index: info.index,
                object: info.object,
            });
            if (event) {
                onFeatureHover(
                    {
                        ...event,
                        coordinateSystem,
                        spatialData,
                        entryId,
                        pickInfo: info,
                    },
                    anchor,
                );
            }
        },
        [coordinateSystem, onFeatureHover, onFeatureHoverEnd, renderer, spatialData],
    );

    if (!stack?.entries.length) {
        return <div className="h-full w-full p-2">Preparing layer stack…</div>;
    }

    if (spatialViewState === null && renderer.hasEnabledLayers) {
        return <div className="h-full w-full p-2">{renderer.isBlocking ? "Loading layer data…" : "Framing view…"}</div>;
    }

    if (!renderer.hasRenderableInputs) {
        return <div className="h-full w-full p-2">No layers to display</div>;
    }

    if (densityGrid && spatialViewState) {
        return (
            <SpatialDensityGrid
                host={densityGrid}
                renderer={renderer}
                renderStack={stackStructure}
                viewState={spatialViewState}
                onViewStateChange={onSpatialViewStateChange}
                onHover={handleHover}
            />
        );
    }

    return (
        <CrosshairArea enabled={isCrosshairShown(config)}>
            <SpatialViewer
                width={width}
                height={height}
                viewState={spatialViewState}
                onViewStateChange={onSpatialViewStateChange}
                layers={renderer.deckLayers}
                layerOrder={renderer.layerOrder}
                renderStack={stackStructure}
                vivLayerProps={renderer.vivLayerProps.length > 0 ? renderer.vivLayerProps : undefined}
                // A prop, not part of `deckProps`: `SpatialViewer` spreads `deckProps` and then
                // sets `onHover` from this prop, so one passed inside `deckProps` is overwritten
                // with `undefined` and feature tooltips never appear. Only `info`: deck also
                // passes its event, which `handleHover` would take for the tooltip's anchor.
                onHover={(info) => handleHover(info)}
                deckProps={deckProps}
                deckRef={deckRef}
            />
        </CrosshairArea>
    );
});

/** The chart's region grid, with the render stack read as the single view reads it. */
const SpatialRegionGridFromRenderStack = observer(function SpatialRegionGridFromRenderStack({
    host,
    layerInputs,
    pointsLayer,
}: {
    host: SpatialRegionGridHost;
    layerInputs: RenderStackLayerInputs;
    pointsLayer: Layer | undefined;
}) {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const stack = config.renderStack;
    void chart.renderStackGeneration;
    const stackStructure = useRenderStackStructure(stack);
    const vivPassthrough = useSpatialVivPassthrough(stack);
    const pointsEntry = stack && findRenderStackEntry(stack, deckHostLayerId("scatter"));
    return (
        <SpatialRegionGrid
            host={host}
            layerInputs={layerInputs}
            renderStack={stackStructure}
            vivPassthrough={vivPassthrough}
            pointsLayer={pointsLayer}
            pointsVisible={pointsEntry?.visible !== false}
        />
    );
});

function SpatialDataChartRoot() {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const { vivStores } = chart;
    return (
        <VivProvider vivStores={vivStores}>
            <SpatialDataMainChart />
        </VivProvider>
    );
}

const SpatialDataMainChart = observer(() => {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const [hoveredContour, setHoveredContour] = useState<HoveredContour | null>(null);
    useEffect(() => {
        ensureParquetWorker();
    }, []);

    return (
        <SpatialAnnotationProvider
            chart={chart}
            hoveredContour={hoveredContour}
            setHoveredContour={setHoveredContour}
            separateDensityLayer
        >
            <SpatialDataViewer setHoveredContour={setHoveredContour} />
        </SpatialAnnotationProvider>
    );
});

/**
 * What the canvas shows when it throws.
 *
 * `FallbackComponent` is handed `FallbackProps` — `{ error, resetErrorBoundary }` —
 * not the error itself, and `error` is `unknown` because anything can be thrown.
 * Both matter: the details dialog renders `error.message` as a React child, so an
 * `Error` object reaching that field makes OPENING the dialog throw its own
 * "Objects are not valid as a React child", replacing the report with a second
 * failure. Narrow to strings here, and pass the stack — it is the only part worth
 * having, and nothing was passing it.
 */
function SpatialCanvasErrorFallback({ error, layers }: FallbackProps & { layers: RenderStackLayerInputs["layers"] }) {
    const isError = error instanceof Error;
    return (
        <ErrorComponentReactWrapper
            error={{
                message: isError ? error.message : String(error),
                ...(isError && error.stack ? { stack: error.stack } : {}),
            }}
            extraMetaData={summariseLayersForError(layers)}
            title="Error showing spatialdata canvas"
        />
    );
}

/**
 * Which layers were being drawn, small enough to put on a clipboard.
 *
 * The dialog `JSON.stringify`s this. A layer config is not safe to stringify
 * whole: `featureState.fillColorByFeatureId` carries one entry per feature, so on
 * a real dataset "copy error details" would serialise hundreds of thousands of
 * keys and hang the tab — turning a report into a second failure. What actually
 * helps diagnosis is which element each layer pointed at and what it was being
 * coloured by.
 */
function summariseLayersForError(layers: RenderStackLayerInputs["layers"]) {
    return {
        layers: Object.entries(layers).map(([id, layer]) => ({
            id,
            type: layer?.type,
            elementKey: layer?.elementKey,
            ...(layer && "fillColorByColumn" in layer && layer.fillColorByColumn
                ? {
                      fillColorByColumn: layer.fillColorByColumn.columnName,
                      fillColorMode: layer.fillColorByColumn.mode,
                  }
                : {}),
        })),
    };
}

const SpatialDataViewer = observer(
    ({ setHoveredContour }: { setHoveredContour: (contour: HoveredContour | null) => void }) => {
        const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
        const config = useConfig<SpatialDataMdvReactConfig>();
        const rawRegion = useRegion();
        const region = getSpatialRegionMetadata(rawRegion);
        const coordinateSystem = region?.spatial?.coordinate_system ?? null;
        const { spatialData, loading, error } = useChartSpatialData();
        const viewerStore = useViewerStoreApi();
        const viewState = useViewerStore((store) => store.viewState);
        const spatialViewState = toSpatialViewState(viewState);
        const [width, height] = useChartSize();
        const id = useChartID();
        const deckContainerRef = useRef<HTMLDivElement | null>(null);
        const { scatterProps, selectionLayer } = useSpatialLayers();
        const { scatterplotLayer, greyScatterplotLayer, densityLayer, getTooltip, setScatterKeyboardActive } =
            scatterProps;
        const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();
        const contourConfig = useConfig<DualContourLegacyConfig>();
        const legendFields = useFieldContourLegend(contourConfig.densityFields);
        const legendCategories = useCategoryContourLegend(contourConfig);
        const showLegend = contourConfig.field_legend.display;
        const legendPosition = { x: 10, y: 10 };
        const { device: deckDevice, fallBack: fallBackToWebGL } = useDeckDevice(
            resolveDeckBackend(config.render_backend),
        );
        const activeBackend = deckDevice.status === "ready" ? deckDevice.backend : null;
        const fallbackReason = deckDevice.status === "ready" ? deckDevice.fallbackReason : undefined;
        useEffect(() => {
            chart.setActiveBackend(activeBackend, fallbackReason);
            return () => chart.setActiveBackend(null);
        }, [chart, activeBackend, fallbackReason]);
        // Layer errors (deck's `onError`) and WebGPU device errors (luma's) otherwise only
        // reach the console, and a device error repeats every frame. The tracker ties them
        // to render-stack entries, which the chart turns off and reports.
        //
        // On WebGPU, a chart following the default backend instead moves to WebGL with the
        // layer still on: the failure is most likely WebGPU's, and a failure on WebGL is
        // reported as usual. A chart set to WebGPU explicitly stays on it.
        const failureTracker = useMemo(
            () =>
                new LayerFailureTracker({
                    entryIds: () => (chart.config.renderStack ? getRenderStackEntryIds(chart.config.renderStack) : []),
                    fail: (entryId, failure) => {
                        const onWebGPU = failure.source === "webgpu" || chart.activeBackend === "webgpu";
                        if (onWebGPU && chart.config.render_backend == null) {
                            // after the switch, the old WebGPU device's errors are stale
                            if (chart.activeBackend !== "webgpu") return;
                            const entry =
                                chart.config.renderStack && findRenderStackEntry(chart.config.renderStack, entryId);
                            const name = entry ? renderStackEntryDisplayName(entry, chart.config.renderStack) : entryId;
                            console.error(`[MDV] ${name} failed on WebGPU; using WebGL:`, failure.reason);
                            fallBackToWebGL(`${name} failed on WebGPU: ${failure.reason}`);
                            return;
                        }
                        chart.recordLayerFailure(entryId, failure);
                    },
                    canvas: (message, stack) => chart.recordCanvasError(message, stack),
                }),
            [chart, fallBackToWebGL],
        );
        useEffect(() => () => failureTracker.dispose(), [failureTracker]);
        const deckErrorHandlers = useMemo<DeckErrorHandlers>(
            () => ({
                device: (error) => failureTracker.device(error),
                layer: (error, layerId) => failureTracker.layer(error, layerId),
                initialized: (device) => chart.setGpuDescription(describeGpu(device.info)),
            }),
            [chart, failureTracker],
        );
        const {
            deviceProps,
            onDeviceInitialized,
            onError: onDeckError,
        } = useDeckDeviceProps(deckDevice, fallBackToWebGL, deckErrorHandlers);

        useViewStateLink();

        useEffect(() => {
            const saved = config.viewState;
            if (!saved?.target || viewState) return;
            viewerStore.setState({
                viewState: {
                    target: [...saved.target],
                    zoom: saved.zoom,
                } as OrthographicViewState,
            });
        }, [config.viewState, viewState, viewerStore]);

        useEffect(() => {
            if (scatterProps.viewState) {
                viewerStore.setState({ viewState: scatterProps.viewState });
            }
        }, [scatterProps.viewState, viewerStore]);

        useEffect(() => {
            if (loading || !spatialData || !coordinateSystem) return;
            seedRenderStackFromSpatialData(config, spatialData, coordinateSystem, chart);
        }, [chart, config, coordinateSystem, loading, spatialData]);

        const storeUrl = spatialData?.url ?? null;
        useEffect(() => {
            chart.setSpatialDataUrl(storeUrl);
            return () => chart.setSpatialDataUrl(null);
        }, [chart, storeUrl]);

        useProjectMdvFieldSpecs(config, chart);

        const deckOverlaySources = useMemo<MdvDeckOverlayLayers>(
            () => ({
                grey_scatter: greyScatterplotLayer,
                scatter: scatterplotLayer,
                density: densityLayer,
                gate_display: gateDisplayLayer,
                selection: selectionLayer,
                gate_labels: gateLabelLayer,
            }),
            [gateLabelLayer, gateDisplayLayer, scatterplotLayer, greyScatterplotLayer, densityLayer, selectionLayer],
        );

        const hostLayerResolver = useMemo(() => createMdvHostLayerResolver(deckOverlaySources), [deckOverlaySources]);

        const {
            layers: baseLayers,
            layerOrder: baseLayerOrder,
            deckLayers,
        } = useRenderStackAdapter({
            stack: config.renderStack,
            generation: chart.renderStackGeneration,
            propsGeneration: chart.renderStackPropsGeneration,
            hostLayerResolver,
        });
        const { layers, layerOrder, mdvTooltipItems } = useAssociatedShapesLayerInputs(
            spatialData ?? undefined,
            { layers: baseLayers, layerOrder: baseLayerOrder },
            region?.spatial?.file,
        );

        const onSpatialViewStateChange = useCallback(
            (next: SpatialCanvasViewState) => {
                viewerStore.setState({
                    viewState: toMdvViewState(next, viewState),
                });
            },
            [viewerStore, viewState],
        );

        const getTooltipContent = useCallback(
            (info: PickingInfo) =>
                getCombinedScatterTooltip(info, {
                    // the render stack draws the gate layers under its own ids, not the layers'
                    gateDisplayLayerId: gateDisplayLayer ? deckHostLayerId("gate_display") : undefined,
                    gateLabelLayerId: gateLabelLayer ? deckHostLayerId("gate_labels") : undefined,
                    getPointTooltip: getTooltip,
                }),
            [gateDisplayLayer, gateLabelLayer, getTooltip],
        );

        const {
            clearTooltip,
            getTooltip: getPortalTooltip,
            suppressTooltipUntilPointerUp,
            tooltipPortal,
        } = useOuterContainerDeckTooltip(getTooltipContent, deckContainerRef);
        const outerContainer = useOuterContainer();
        const [featureTooltip, setFeatureTooltip] = useState<{
            html: string;
            x: number;
            y: number;
        } | null>(null);

        const onFeatureHover = useCallback(
            (event: SpatialFeaturePickEvent, anchorElement?: HTMLElement) => {
                const info = event.pickInfo;
                // The viewer fills in the fields it can read from the table's obs; MDV adds
                // the rest from the datasource row (see `splitTooltipFields`).
                const extra =
                    event.featureId === undefined ? [] : mdvTooltipItems(event.entryId, String(event.featureId));
                const tooltip =
                    extra.length === 0
                        ? event.tooltip
                        : {
                              ...(event.tooltip ?? { title: String(event.featureId) }),
                              items: [...(event.tooltip?.items ?? []), ...extra],
                          };
                if (!tooltip || !Number.isFinite(info.x) || !Number.isFinite(info.y)) {
                    setFeatureTooltip(null);
                    return;
                }
                const anchor = anchorElement ?? deckContainerRef.current;
                if (!anchor) return;
                const rect = anchor.getBoundingClientRect();
                setFeatureTooltip({
                    html: formatSpatialFeatureTooltipHtml(tooltip),
                    x: rect.left + info.x,
                    y: rect.top + info.y,
                });
            },
            [mdvTooltipItems],
        );

        /** Deck reports moving off a feature as an unpicked hover, not as an event. */
        const onFeatureHoverEnd = useCallback(() => setFeatureTooltip(null), []);

        const featureTooltipPortal =
            featureTooltip && outerContainer
                ? createPortal(
                      <div
                          style={{
                              position: "fixed",
                              left: featureTooltip.x + 12,
                              top: featureTooltip.y + 12,
                              zIndex: 10000,
                              pointerEvents: "none",
                              maxWidth: 300,
                              borderRadius: 6,
                              border: "1px solid rgba(255,255,255,0.1)",
                              backgroundColor: "rgba(10, 10, 10, 0.9)",
                              padding: "8px 10px",
                              fontSize: 12,
                              color: "#f5f5f5",
                          }}
                          // biome-ignore lint/security/noDangerouslySetInnerHtml: html is sanitized with DOMPurify in formatSpatialFeatureTooltipHtml
                          dangerouslySetInnerHTML={{ __html: featureTooltip.html }}
                      />,
                      outerContainer,
                  )
                : null;

        const deckProps: Partial<DeckGLProps> = useMemo(
            () => ({
                getTooltip: getPortalTooltip,
                id: `${id}spatialdata-deck`,
                controller: {
                    type: GatedOrthographicController,
                    doubleClickZoom: false,
                    ...controllerOptions,
                },
                ...(deviceProps && { deviceProps }),
                onDeviceInitialized,
                onError: onDeckError,
            }),
            [id, getPortalTooltip, controllerOptions, deviceProps, onDeviceInitialized, onDeckError],
        );

        const regionGridShown = isRegionGridShown(config);
        const regionGridLayerInputs = useMemo(
            () => ({ layers: baseLayers, layerOrder: baseLayerOrder }),
            [baseLayers, baseLayerOrder],
        );
        const regionGrid = useMemo<SpatialRegionGridHost | undefined>(
            () =>
                regionGridShown
                    ? {
                          device: deckDevice,
                          deviceProps,
                          onDeviceInitialized,
                          onError: onDeckError,
                          controllerOptions,
                          getTooltipContent,
                          setKeyboardActive: setScatterKeyboardActive,
                          onPointerLeave: onFeatureHoverEnd,
                      }
                    : undefined,
            [
                regionGridShown,
                deckDevice,
                deviceProps,
                onDeviceInitialized,
                onDeckError,
                controllerOptions,
                getTooltipContent,
                setScatterKeyboardActive,
                onFeatureHoverEnd,
            ],
        );

        const densityGridCells = useDensityGridCells();
        // with no field to make a cell of, the one canvas, without an overlay (as the scatter plot)
        const gridShown = isDensityGridShown(config) && densityGridCells.cells.length > 0;
        const densityGrid = useMemo<SpatialDensityGridHost | undefined>(
            () =>
                gridShown
                    ? {
                          ...densityGridCells,
                          device: deckDevice,
                          deviceProps,
                          onDeviceInitialized,
                          onError: onDeckError,
                          controllerOptions,
                          getTooltipContent,
                          setKeyboardActive: setScatterKeyboardActive,
                          onPointerLeave: onFeatureHoverEnd,
                      }
                    : undefined,
            [
                gridShown,
                densityGridCells,
                deckDevice,
                deviceProps,
                onDeviceInitialized,
                onDeckError,
                controllerOptions,
                getTooltipContent,
                setScatterKeyboardActive,
                onFeatureHoverEnd,
            ],
        );

        if (!coordinateSystem || !region?.spatial?.file) {
            return (
                <div className="h-full w-full p-2">
                    SpatialData.js viewer requires region.spatial.file and region.spatial.coordinate_system metadata.
                </div>
            );
        }

        if (error) {
            return (
                <div className="h-full w-full p-2">
                    Failed to load SpatialData store: {error.message}
                    <Button size="small" onClick={() => chart.spatialDataReader.retry()}>
                        Retry
                    </Button>
                </div>
            );
        }

        if (loading || !spatialData) {
            return <div className="h-full w-full p-2">Loading SpatialData store…</div>;
        }

        if (!config.renderStack?.entries.length) {
            return <div className="h-full w-full p-2">Preparing layer stack…</div>;
        }

        if (deckDevice.status === "pending") {
            return <div className="h-full w-full p-2">Starting WebGPU…</div>;
        }

        return (
            <Profiler id="spatial.viewer" onRender={onSpatialProfilerRender}>
                {/* selection tools, gates and density are not in region cells yet */}
                {!regionGrid && <SelectionOverlay />}
                {/* in the grid each cell is labelled with its field; categories still need the legend */}
                {showLegend && !regionGrid && (!densityGrid || legendCategories.items.length > 0) && (
                    <FieldContourLegend
                        fields={densityGrid ? [] : legendFields}
                        categories={legendCategories.items}
                        categoryLabel={legendCategories.label}
                        position={legendPosition}
                        onContourHover={setHoveredContour}
                    />
                )}
                <div
                    ref={deckContainerRef}
                    aria-label="SpatialData.js image viewer"
                    style={{ width: "100%", height: "100%", outline: "none" }}
                    onPointerDown={suppressTooltipUntilPointerUp}
                    onMouseDown={() => {
                        setScatterKeyboardActive(true);
                    }}
                    onMouseEnter={() => setScatterKeyboardActive(true)}
                    onMouseLeave={() => {
                        clearTooltip();
                        setFeatureTooltip(null);
                        setScatterKeyboardActive(false);
                    }}
                >
                    <div style={{ width, height, position: "relative" }}>
                        <Profiler id="spatial.canvas" onRender={onSpatialProfilerRender}>
                            <ErrorBoundary
                                // deck picks its device once, when created: remount on a backend change.
                                key={deckDevice.backend}
                                // On WebGPU, a chart following the default retries on WebGL before
                                // showing the error; a failure there too is still reported.
                                onError={(error) => {
                                    if (deckDevice.status !== "ready" || deckDevice.backend !== "webgpu") return;
                                    if (config.render_backend != null) return;
                                    const message = error instanceof Error ? error.message : String(error);
                                    fallBackToWebGL(`the canvas failed on WebGPU: ${message}`);
                                }}
                                fallbackRender={(props) => <SpatialCanvasErrorFallback {...props} layers={layers} />}
                                // Retry when the layer set actually changes, so a canvas
                                // felled by one bad layer recovers once it is removed or
                                // reconfigured instead of staying dead until remount.
                                // `layers` keeps its identity across cosmetic edits (see
                                // `layerConfigReplacementSignature`), so an opacity drag
                                // cannot spin this into a re-throw loop.
                                resetKeys={[layers, layerOrder]}
                            >
                                {regionGrid ? (
                                    <SpatialRegionGridFromRenderStack
                                        host={regionGrid}
                                        layerInputs={regionGridLayerInputs}
                                        pointsLayer={scatterplotLayer}
                                    />
                                ) : (
                                    <SpatialCanvasFromRenderStack
                                        spatialData={spatialData}
                                        coordinateSystem={coordinateSystem}
                                        spatialViewState={spatialViewState}
                                        onSpatialViewStateChange={onSpatialViewStateChange}
                                        deckLayers={deckLayers}
                                        layers={layers}
                                        layerOrder={layerOrder}
                                        width={width}
                                        height={height}
                                        deckProps={deckProps}
                                        onFeatureHover={onFeatureHover}
                                        onFeatureHoverEnd={onFeatureHoverEnd}
                                        densityGrid={densityGrid}
                                    />
                                )}
                            </ErrorBoundary>
                        </Profiler>
                    </div>
                </div>
                {/* the grid has its own, over its scrolling cells */}
                {!densityGrid && !regionGrid && <ViewGestureHint />}
                {tooltipPortal}
                {featureTooltipPortal}
            </Profiler>
        );
    },
);

export default SpatialDataChartRoot;
