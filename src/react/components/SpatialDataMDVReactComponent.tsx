import {
    SpatialCanvasProvider,
    type ViewState as SpatialCanvasViewState,
} from "@spatialdata/vis";
import { SpatialDataProvider, useSpatialData } from "@spatialdata/react";
import type { Layer } from "@deck.gl/core";
import type {
    DeckGLProps,
    OrbitViewState,
    OrthographicViewState,
    PickingInfo,
} from "deck.gl";
import { observer } from "mobx-react-lite";
import { runInAction } from "mobx";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";

import { getProjectURL } from "@/dataloaders/DataLoaderUtil";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import type { FieldName } from "@/charts/charts";
import {
    VivProvider,
    useViewerStore,
    useViewerStoreApi,
} from "./avivatorish/state";
import SelectionOverlay from "./SelectionOverlay";
import FieldContourLegend from "./FieldContourLegend";
import { useFieldContourLegend } from "../contour_state";
import type { DualContourLegacyConfig } from "../contour_state";
import { useChart } from "../context";
import {
    useChartID,
    useChartSize,
    useConfig,
    useRegion,
} from "../hooks";
import { useViewStateLink } from "../chartLinkHooks";
import {
    SpatialAnnotationProvider,
    useSpatialLayers,
} from "../spatial_context";
import useGateLayers from "../hooks/useGateLayers";
import { useOuterContainerDeckTooltip } from "../hooks/useOuterContainerDeckTooltip";
import {
    type DeckOverlayId,
    type SpatialLayerStackConfig,
    applySpatialLayerStack,
    canvasLayerOrder,
    createDefaultSpatialLayerStack,
    deckStackKey,
    isDeckStackKey,
    normalizeSpatialLayerStack,
} from "@/react/spatial_layer_stack";
import { SpatialCanvasDeckView } from "./SpatialCanvasDeckView";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "./SpatialDataMDVReact";

type SpatialRegionMetadata = {
    spatial?: {
        file?: string;
        coordinate_system?: string;
    };
};

type SpatialCanvas2DViewState = {
    target: [number, number];
    zoom: number;
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

function toSpatialViewState(
    viewState: OrthographicViewState | OrbitViewState | null | undefined,
): SpatialCanvas2DViewState | null {
    if (!viewState || !Array.isArray(viewState.target)) return null;
    const [x, y] = viewState.target;
    if (typeof x !== "number" || typeof y !== "number") return null;
    const zoom = typeof viewState.zoom === "number" ? viewState.zoom : 0;
    return {
        target: [x, y],
        zoom,
    };
}

function toMdvViewState(
    viewState: SpatialCanvasViewState,
    previousViewState: OrthographicViewState | OrbitViewState | null | undefined,
): OrthographicViewState {
    const [x, y] = viewState.target;
    const previousTarget = Array.isArray(previousViewState?.target)
        ? previousViewState.target
        : [0, 0, 0];
    return {
        ...(previousViewState ?? {}),
        target: [x, y, previousTarget[2] ?? 0],
        zoom: viewState.zoom,
    };
}

function withStableDeckOverlayId(layer: Layer, deckId: DeckOverlayId): Layer {
    return layer.clone({ id: deckStackKey(deckId) });
}

function buildDeckOverlayLayers(
    stack: SpatialLayerStackConfig,
    overlays: Record<DeckOverlayId, Layer | null>,
): Layer[] {
    const layers: Layer[] = [];
    for (const key of stack.stackOrder) {
        if (!isDeckStackKey(key)) continue;
        const entry = stack.entries[key];
        if (!entry || entry.kind !== "deck" || !entry.visible) continue;
        const source = overlays[entry.deckId];
        if (!source) continue;
        layers.push(withStableDeckOverlayId(source, entry.deckId));
    }
    return layers;
}

function SpatialDataChartRoot() {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const { vivStores } = chart;
    return (
        <VivProvider vivStores={vivStores}>
            <SpatialCanvasProvider store={chart.spatialCanvasStore}>
                <SpatialDataMainChart />
            </SpatialCanvasProvider>
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

    return (
        <SpatialDataProvider source={spatialDataUrl ?? undefined}>
            <SpatialAnnotationProvider
                chart={chart}
                hoveredFieldId={hoveredField}
                setHoveredFieldId={setHoveredField}
            >
                <SpatialDataViewer
                    setHoveredField={setHoveredField}
                />
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
        const spatialdataPath = region?.spatial?.file ?? null;
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
        const {
            gateLabelLayer,
            gateDisplayLayer,
            controllerOptions,
        } = useGateLayers();
        const contourConfig = useConfig<DualContourLegacyConfig>();
        const legendFields = useFieldContourLegend(contourConfig.densityFields);
        const showLegend = contourConfig.field_legend.display;
        const legendPosition = { x: 10, y: 10 };

        useViewStateLink();

        useEffect(() => {
            if (scatterProps.viewState) {
                viewerStore.setState({ viewState: scatterProps.viewState });
            }
        }, [scatterProps.viewState, viewerStore.setState]);

        const stackSeededRef = useRef(false);
        useEffect(() => {
            if (loading || !spatialData || !coordinateSystem || stackSeededRef.current) return;
            stackSeededRef.current = true;
            runInAction(() => {
                const stack = config.spatialLayerStack;
                if (!stack) return;
                const next =
                    stack.stackOrder.length === 0
                        ? createDefaultSpatialLayerStack(spatialData, coordinateSystem)
                        : normalizeSpatialLayerStack(stack, spatialData, coordinateSystem);
                applySpatialLayerStack(stack, next);
                chart.updateCanvasLayerState(stack);
            });
        }, [chart, config, coordinateSystem, loading, spatialData]);

        const stack = config.spatialLayerStack;
        const unifiedLayerOrder = useMemo(
            () => (stack ? canvasLayerOrder(stack) : []),
            [stack],
        );

        const deckOverlaySources = useMemo(
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

        const deckLayers = useMemo(
            () => (stack ? buildDeckOverlayLayers(stack, deckOverlaySources) : []),
            [deckOverlaySources, stack],
        );

        const getTooltipContent = useCallback(
            (info: PickingInfo) => {
                return getCombinedScatterTooltip(
                    info,
                    {
                        gateDisplayLayerId: gateDisplayLayer?.id,
                        gateLabelLayerId: gateLabelLayer?.id,
                        getPointTooltip: getTooltip,
                    },
                );
            },
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

        return (
            <>
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
                        setScatterKeyboardActive(false);
                    }}
                >
                    <div style={{ width, height, position: "relative" }}>
                        <SpatialCanvasDeckView
                            spatialData={spatialData}
                            coordinateSystem={coordinateSystem}
                            spatialdataPath={spatialdataPath}
                            spatialViewState={spatialViewState}
                            width={width}
                            height={height}
                            deckLayers={deckLayers}
                            unifiedLayerOrder={unifiedLayerOrder}
                            deckProps={deckProps}
                            toMdvViewState={toMdvViewState}
                        />
                    </div>
                </div>
                <SpatialDataLink />
                {tooltipPortal}
            </>
        );
    },
);

export default SpatialDataChartRoot;
