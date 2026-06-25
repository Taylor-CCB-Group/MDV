import { ColorPaletteExtension } from "@hms-dbmi/viv";
import { SpatialDataProvider, useSpatialData } from "@spatialdata/react";
import {
    SpatialViewer,
    useSpatialCanvasRendererFromLayerInputs,
    type SpatialFeaturePickEvent,
    type ViewState as SpatialCanvasViewState,
    type VivImageLayerContext,
    type VivImagePropsResolver,
} from "@spatialdata/vis";
import type { DeckGLProps, OrthographicViewState, PickingInfo } from "deck.gl";
import { observer } from "mobx-react-lite";
import { Profiler, useCallback, useEffect, useMemo, useRef, useState } from "react";
import { createPortal } from "react-dom";

import { getProjectURL } from "@/dataloaders/DataLoaderUtil";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import type { FieldName } from "@/charts/charts";
import { createImageLayerRegistry } from "@/react/spatialdata/image_layer_registry";
import { onSpatialProfilerRender } from "@/react/spatialdata/perf";
import {
    createMdvHostLayerResolver,
    useRenderStackAdapter,
    type MdvDeckOverlayLayers,
} from "@/react/spatialdata/render_stack_adapter";
import { seedRenderStackFromSpatialData } from "@/react/spatialdata/render_stack_control";
import { toMdvViewState, toSpatialViewState } from "@/react/spatialdata/view_state_bridge";
import { ensureChunkWorker } from "@/react/spatialdata/ensureChunkWorker";
import { formatSpatialFeatureTooltipHtml } from "@/react/spatialdata/spatial_feature_tooltip";
import VivContrastExtension from "@/webgl/VivContrastExtension";
import { useOuterContainer } from "../screen_state";
import { useViewStateLink } from "../chartLinkHooks";
import { useChart } from "../context";
import { useChartID, useChartSize, useConfig, useRegion } from "../hooks";
import useGateLayers from "../hooks/useGateLayers";
import { useOuterContainerDeckTooltip } from "../hooks/useOuterContainerDeckTooltip";
import { SpatialAnnotationProvider, useSpatialLayers } from "../spatial_context";
import { useFieldContourLegend } from "../contour_state";
import type { DualContourLegacyConfig } from "../contour_state";
import {
    VivProvider,
    useViewerStore,
    useViewerStoreApi,
} from "./avivatorish/state";
import FieldContourLegend from "./FieldContourLegend";
import SelectionOverlay from "./SelectionOverlay";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "./SpatialDataMDVReact";

type SpatialRegionMetadata = {
    spatial?: {
        file?: string;
        coordinate_system?: string;
    };
};

function getSpatialRegionMetadata(region: unknown): SpatialRegionMetadata | null {
    if (!region || typeof region !== "object") return null;
    return region;
}

function getSpatialDataUrl(region: SpatialRegionMetadata) {
    const file = region.spatial?.file;
    if (!file) return null;
    return getProjectURL(`spatial/${file}`);
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
}: {
    spatialData: NonNullable<ReturnType<typeof useSpatialData>["spatialData"]>;
    coordinateSystem: string;
    spatialViewState: ReturnType<typeof toSpatialViewState>;
    onSpatialViewStateChange: (next: SpatialCanvasViewState) => void;
    deckLayers: ReturnType<typeof useRenderStackAdapter>["deckLayers"];
    layers: ReturnType<typeof useRenderStackAdapter>["layers"];
    layerOrder: ReturnType<typeof useRenderStackAdapter>["layerOrder"];
    width: number;
    height: number;
    deckProps: Partial<DeckGLProps>;
    onFeatureHover: (event: SpatialFeaturePickEvent) => void;
}) {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const stack = config.renderStack;
    void chart.renderStackGeneration;

    const extensions = useMemo(
        () => [new ColorPaletteExtension(), new VivContrastExtension()],
        [],
    );

    const vivImagePropsResolver = useCallback<VivImagePropsResolver>(
        (ctx: VivImageLayerContext) => {
            if (!stack) return undefined;
            const entry = stack.entries.find((item) => item.id === ctx.layerId);
            if (!entry || entry.kind !== "spatial") return undefined;
            const saved = entry.props.vivLayerProps;
            if (!saved || typeof saved !== "object") return undefined;
            return saved as Record<string, unknown>;
        },
        [stack],
    );

    const vivPassthrough = useMemo(
        () => ({
            vivImageExtensions: extensions,
            vivImageExtensionResolver: () => extensions,
            vivImagePropsResolver,
        }),
        [extensions, vivImagePropsResolver],
    );

    const renderer = useSpatialCanvasRendererFromLayerInputs({
        spatialData,
        coordinateSystem,
        layerInputs: { layers, layerOrder },
        renderOrder: layerOrder,
        viewState: spatialViewState,
        onViewStateChange: onSpatialViewStateChange,
        width,
        height,
        hostDeckLayers,
        sortDeckLayers: true,
        autoFit: spatialViewState === null,
        vivPassthrough,
    });

    useEffect(() => {
        if (!stack) {
            chart.setImageLayerRegistry(undefined);
            return;
        }
        chart.setImageLayerRegistry(
            createImageLayerRegistry(
                stack,
                renderer.getImageLoadedDataByElementKey,
                renderer.getLayerLoadState,
            ),
        );
        return () => chart.setImageLayerRegistry(undefined);
    }, [
        chart,
        stack,
        chart.renderStackGeneration,
        renderer.getImageLoadedDataByElementKey,
        renderer.getLayerLoadState,
        renderer.isBlocking,
        renderer.isLoading,
    ]);

    const handleHover = useCallback(
        (info: PickingInfo) => {
            if (!info.picked || typeof info.x !== "number" || typeof info.y !== "number") {
                return;
            }
            const layerId = (typeof info.layer?.id === "string" ? info.layer.id : "").replace(
                /-#.*#$/,
                "",
            );
            const event = renderer.getFeaturePickEvent(layerId, {
                index: info.index,
                object: info.object,
            });
            if (event) {
                onFeatureHover({
                    ...event,
                    coordinateSystem,
                    spatialData,
                    pickInfo: info,
                });
            }
        },
        [coordinateSystem, onFeatureHover, renderer, spatialData],
    );

    const mergedDeckProps = useMemo(
        () => ({
            ...deckProps,
            onHover: handleHover,
        }),
        [deckProps, handleHover],
    );

    if (!stack?.entries.length) {
        return <div className="h-full w-full p-2">Preparing layer stack…</div>;
    }

    if (spatialViewState === null && renderer.hasEnabledLayers) {
        return (
            <div className="h-full w-full p-2">
                {renderer.isBlocking ? "Loading layer data…" : "Framing view…"}
            </div>
        );
    }

    if (!renderer.hasRenderableInputs) {
        return <div className="h-full w-full p-2">No layers to display</div>;
    }

    return (
        <SpatialViewer
            width={width}
            height={height}
            viewState={spatialViewState}
            onViewStateChange={onSpatialViewStateChange}
            layers={renderer.deckLayers}
            layerOrder={renderer.layerOrder}
            vivLayerProps={renderer.vivLayerProps.length > 0 ? renderer.vivLayerProps : undefined}
            deckProps={mergedDeckProps}
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

function SpatialDataLink() {
    const { spatialData } = useSpatialData();
    const spatialDataUrl = spatialData?.url;
    if (!spatialDataUrl) return null;
    const demoUrl = `https://taylor-ccb-group.github.io/SpatialData.js/docs/demo/?url=${encodeURIComponent(spatialDataUrl)}`;
    return (
        <div className="legend-container pointer-events-auto absolute bottom-2 right-2 z-[4] rounded-tl border-[0.5px] border-current">
            <a
                href={demoUrl}
                target="_blank"
                rel="noopener noreferrer"
                className="text-xs text-current no-underline hover:underline whitespace-nowrap"
            >
                Open in spatialdata.js
            </a>
        </div>
    );
}

const SpatialDataMainChart = observer(() => {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const [hoveredField, setHoveredField] = useState<FieldName | null>(null);
    const rawRegion = useRegion();
    const region = getSpatialRegionMetadata(rawRegion);
    const spatialDataUrl = region ? getSpatialDataUrl(region) : null;
    ensureChunkWorker();

    return (
        <SpatialDataProvider source={spatialDataUrl ?? undefined}>
            <SpatialAnnotationProvider
                chart={chart}
                hoveredFieldId={hoveredField}
                setHoveredFieldId={setHoveredField}
            >
                <SpatialDataViewer setHoveredField={setHoveredField} />
            </SpatialAnnotationProvider>
        </SpatialDataProvider>
    );
});

const SpatialDataViewer = observer(
    ({
        setHoveredField,
    }: {
        setHoveredField: (fieldId: FieldName | null) => void;
    }) => {
        const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
        const config = useConfig<SpatialDataMdvReactConfig>();
        const rawRegion = useRegion();
        const region = getSpatialRegionMetadata(rawRegion);
        const coordinateSystem = region?.spatial?.coordinate_system ?? null;
        const { spatialData, loading, error } = useSpatialData();
        const viewerStore = useViewerStoreApi();
        const viewState = useViewerStore((store) => store.viewState);
        const spatialViewState = toSpatialViewState(viewState);
        const [width, height] = useChartSize();
        const id = useChartID();
        const deckContainerRef = useRef<HTMLDivElement | null>(null);
        const { scatterProps, selectionLayer } = useSpatialLayers();
        const {
            scatterplotLayer,
            greyScatterplotLayer,
            getTooltip,
            setScatterKeyboardActive,
        } = scatterProps;
        const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();
        const contourConfig = useConfig<DualContourLegacyConfig>();
        const legendFields = useFieldContourLegend(contourConfig.densityFields);
        const showLegend = contourConfig.field_legend.display;
        const legendPosition = { x: 10, y: 10 };

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

        const deckOverlaySources = useMemo<MdvDeckOverlayLayers>(
            () => ({
                grey_scatter: greyScatterplotLayer,
                scatter: scatterplotLayer,
                gate_display: gateDisplayLayer,
                selection: selectionLayer,
                gate_labels: gateLabelLayer,
            }),
            [
                gateLabelLayer,
                gateDisplayLayer,
                scatterplotLayer,
                greyScatterplotLayer,
                selectionLayer,
            ],
        );

        const hostLayerResolver = useMemo(
            () => createMdvHostLayerResolver(deckOverlaySources),
            [deckOverlaySources],
        );

        const { layers, layerOrder, deckLayers } = useRenderStackAdapter({
            stack: config.renderStack,
            generation: chart.renderStackGeneration,
            hostLayerResolver,
        });

        const onSpatialViewStateChange = useCallback(
            (next: SpatialCanvasViewState) => {
                viewerStore.setState({
                    viewState: toMdvViewState(next, viewState),
                });
            },
            [viewerStore, viewState],
        );

        const getTooltipContent = useCallback(
            (info: PickingInfo) => {
                return getCombinedScatterTooltip(info, {
                    gateDisplayLayerId: gateDisplayLayer?.id,
                    gateLabelLayerId: gateLabelLayer?.id,
                    getPointTooltip: getTooltip,
                });
            },
            [gateDisplayLayer?.id, gateLabelLayer?.id, getTooltip],
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
            (event: SpatialFeaturePickEvent) => {
                const info = event.pickInfo;
                if (!event.tooltip || !Number.isFinite(info.x) || !Number.isFinite(info.y)) {
                    setFeatureTooltip(null);
                    return;
                }
                const anchor = deckContainerRef.current;
                if (!anchor) return;
                const rect = anchor.getBoundingClientRect();
                setFeatureTooltip({
                    html: formatSpatialFeatureTooltipHtml(event.tooltip),
                    x: rect.left + info.x,
                    y: rect.top + info.y,
                });
            },
            [],
        );

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
                    doubleClickZoom: false,
                    dragPan: controllerOptions.dragPan,
                },
            }),
            [id, getPortalTooltip, controllerOptions],
        );

        if (!coordinateSystem || !region?.spatial?.file) {
            return (
                <div className="h-full w-full p-2">
                    SpatialData.js viewer requires region.spatial.file and
                    region.spatial.coordinate_system metadata.
                </div>
            );
        }

        if (error) {
            return (
                <div className="h-full w-full p-2">
                    Failed to load SpatialData store: {error.message}
                </div>
            );
        }

        if (loading || !spatialData) {
            return <div className="h-full w-full p-2">Loading SpatialData store…</div>;
        }

        if (!config.renderStack?.entries.length) {
            return <div className="h-full w-full p-2">Preparing layer stack…</div>;
        }

        return (
            <Profiler id="spatial.viewer" onRender={onSpatialProfilerRender}>
                <SelectionOverlay />
                {showLegend && legendFields.length > 0 && (
                    <FieldContourLegend
                        fields={legendFields}
                        position={legendPosition}
                        onFieldHover={setHoveredField}
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
                            />
                        </Profiler>
                    </div>
                </div>
                <SpatialDataLink />
                {tooltipPortal}
                {featureTooltipPortal}
            </Profiler>
        );
    },
);

export default SpatialDataChartRoot;
