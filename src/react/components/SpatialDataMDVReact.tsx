import { readZarr, type SpatialData } from "@spatialdata/core";
import {
    SpatialCanvasProvider,
    SpatialViewer,
    createSpatialCanvasStore,
    useSpatialCanvasRenderer,
    type SpatialCanvasStoreApi,
    type ViewState as SpatialCanvasViewState,
} from "@spatialdata/vis";
import type { Layer } from "@deck.gl/core";
import type {
    DeckGLProps,
    OrbitViewState,
    OrthographicViewState,
    PickingInfo,
} from "deck.gl";
import { observer } from "mobx-react-lite";
import { action, makeObservable, observable, runInAction } from "mobx";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";

import BaseChart from "../../charts/BaseChart";
import type DataStore from "@/datastore/DataStore";
import { getProjectURL } from "@/dataloaders/DataLoaderUtil";
import { allNumeric } from "@/lib/columnTypeHelpers";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import type { FieldName } from "@/charts/charts";
import { BaseReactChart } from "./BaseReactChart";
import "../../charts/VivScatterPlot";
import {
    type VivContextType,
    VivProvider,
    applyDefaultChannelState,
    createVivStores,
    useViewerStore,
    useViewerStoreApi,
} from "./avivatorish/state";
import SelectionOverlay from "./SelectionOverlay";
import FieldContourLegend from "./FieldContourLegend";
import { useFieldContourLegend } from "../contour_state";
import type { DualContourLegacyConfig } from "../contour_state";
import type { VivMdvReactConfig, VivRoiConfig } from "./VivMDVReact";
import { useChart, useDataStore } from "../context";
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
import { getSharedScatterSettings } from "./sharedScatterSettings";
import { useOuterContainer } from "../screen_state";
import { scatterDefaults } from "../scatter_state";
import { g, toArray } from "@/lib/utils";
import {
    type DeckOverlayId,
    type SpatialLayerStackConfig,
    applySpatialLayerStack,
    canvasLayerOrder,
    createDefaultSpatialLayerStack,
    deckStackKey,
    isDeckStackKey,
    normalizeSpatialLayerStack,
    stackToCanvasState,
} from "@/react/spatial_layer_stack";
import { inferTableAssociation } from "@/react/spatial_table_association";
import {
    SpatialCanvasRendererProvider,
    type SpatialCanvasRendererValue,
} from "@/react/spatial_canvas_renderer_context";
import SpatialLayerDialogReactWrapper from "./SpatialLayerDialogReactWrapper";

type SpatialRegionMetadata = {
    spatial?: {
        file?: string;
        coordinate_system?: string;
    };
};

type LoadState =
    | { status: "idle"; spatialData: null; error: null }
    | { status: "loading"; spatialData: null; error: null }
    | { status: "loaded"; spatialData: SpatialData; error: null }
    | { status: "error"; spatialData: null; error: Error };

export type SpatialDataMdvReactConfig = VivMdvReactConfig & {
    spatialLayerStack?: SpatialLayerStackConfig;
};

const spatialDataCache = new Map<string, Promise<SpatialData>>();
type SpatialCanvas2DViewState = {
    target: [number, number];
    zoom: number;
};

function getSpatialData(url: string) {
    const existing = spatialDataCache.get(url);
    if (existing) return existing;
    const promise = readZarr(url);
    spatialDataCache.set(url, promise);
    return promise;
}

function getSpatialRegionMetadata(region: unknown): SpatialRegionMetadata | null {
    if (!region || typeof region !== "object") return null;
    return region;
}

function getSpatialDataUrl(region: SpatialRegionMetadata) {
    const file = region.spatial?.file;
    if (!file) return null;
    return getProjectURL(`spatial/${file}`);
}

function useSpatialData(region: SpatialRegionMetadata | null) {
    const url = region ? getSpatialDataUrl(region) : null;
    const [state, setState] = useState<LoadState>({
        status: "idle",
        spatialData: null,
        error: null,
    });

    useEffect(() => {
        if (!url) {
            setState({ status: "idle", spatialData: null, error: null });
            return;
        }
        let cancelled = false;
        setState({ status: "loading", spatialData: null, error: null });
        getSpatialData(url)
            .then((spatialData) => {
                if (!cancelled) {
                    setState({ status: "loaded", spatialData, error: null });
                }
            })
            .catch((error: unknown) => {
                if (cancelled) return;
                setState({
                    status: "error",
                    spatialData: null,
                    error: error instanceof Error
                        ? error
                        : new Error(String(error)),
                });
            });
        return () => {
            cancelled = true;
        };
    }, [url]);

    return { ...state, url };
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

const SpatialDataMainChart = observer(() => {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const [hoveredField, setHoveredField] = useState<FieldName | null>(null);

    return (
        <SpatialAnnotationProvider
            chart={chart}
            hoveredFieldId={hoveredField}
            setHoveredFieldId={setHoveredField}
        >
            <SpatialDataViewer
                setHoveredField={setHoveredField}
            />
        </SpatialAnnotationProvider>
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
        const dataStore = useDataStore();
        const rawRegion = useRegion();
        const region = getSpatialRegionMetadata(rawRegion);
        const coordinateSystem = region?.spatial?.coordinate_system ?? null;
        const spatialdataPath = region?.spatial?.file ?? null;
        const { spatialData, status, error } = useSpatialData(region);
        const viewerStore = useViewerStoreApi();
        const viewState = useViewerStore((store) => store.viewState);
        const spatialViewState = toSpatialViewState(viewState);
        const [width, height] = useChartSize();
        const id = useChartID();
        const deckContainerRef = useRef<HTMLDivElement | null>(null);
        const outerContainer = useOuterContainer();
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
            if (status !== "loaded" || !coordinateSystem || stackSeededRef.current) return;
            stackSeededRef.current = true;
            runInAction(() => {
                const stack = config.spatialLayerStack;
                if (!stack) return;
                const next =
                    stack.stackOrder.length === 0
                        ? createDefaultSpatialLayerStack(spatialData, coordinateSystem)
                        : normalizeSpatialLayerStack(stack, spatialData, coordinateSystem);
                applySpatialLayerStack(stack, next);
                chart.syncSpatialLayerStack(stack);
            });
        }, [chart, config, coordinateSystem, spatialData, status]);

        const stack = config.spatialLayerStack;
        const { layers, layerOrder } = stack
            ? stackToCanvasState(stack)
            : { layers: {}, layerOrder: [] };
        const unifiedLayerOrder = stack ? canvasLayerOrder(stack) : [];

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

        const deckLayers = stack
            ? buildDeckOverlayLayers(stack, deckOverlaySources)
            : [];

        const renderer = useSpatialCanvasRenderer({
            spatialData,
            coordinateSystem,
            layers,
            layerOrder,
            viewState: spatialViewState,
            onViewStateChange: (nextViewState) => {
                viewerStore.setState({
                    viewState: toMdvViewState(
                        nextViewState,
                        viewerStore.getState().viewState,
                    ),
                });
            },
            width,
            height,
            deckLayers,
            autoFit: true,
        });

        const rendererContextValue = useMemo(
            () => ({
                spatialData,
                coordinateSystem,
                availableElements: renderer.availableElements,
                getImageLayerLoadedData: renderer.getImageLayerLoadedData,
                getLabelsLayerLoadedData: renderer.getLabelsLayerLoadedData,
                getLayerLoadState: renderer.getLayerLoadState,
                isLoading: renderer.isLoading,
                isBlocking: renderer.isBlocking,
                inferTableAssociation: (elementKey?: string | null) =>
                    inferTableAssociation({
                        spatialData,
                        spatialdataPath,
                        regionId: config.region,
                        elementKey,
                        dataStore,
                    }),
            }),
            [
                coordinateSystem,
                config.region,
                dataStore,
                renderer.availableElements,
                renderer.getImageLayerLoadedData,
                renderer.getLabelsLayerLoadedData,
                renderer.getLayerLoadState,
                renderer.isBlocking,
                renderer.isLoading,
                spatialData,
                spatialdataPath,
            ],
        );

        useEffect(() => {
            chart.spatialRendererContext = rendererContextValue;
        }, [chart, rendererContextValue]);

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

        if (status === "error") {
            return (
                <div className="h-full w-full p-2">
                    Failed to load SpatialData store: {error.message}
                </div>
            );
        }

        return (
            <SpatialCanvasRendererProvider value={rendererContextValue}>
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
                        {renderer.isBlocking && (
                            <div className="absolute right-2 top-2 rounded bg-black/70 px-2 py-1 text-[11px] text-white">
                                Loading layer data...
                            </div>
                        )}
                        {!renderer.isBlocking && renderer.isLoading && (
                            <div className="absolute right-2 top-2 rounded bg-neutral-900/80 px-2 py-1 text-[11px] text-white">
                                Refreshing layer metadata...
                            </div>
                        )}
                        <SpatialViewer
                            width={width}
                            height={height}
                            viewState={spatialViewState}
                            onViewStateChange={(nextViewState) => {
                                viewerStore.setState({
                                    viewState: toMdvViewState(
                                        nextViewState,
                                        viewerStore.getState().viewState,
                                    ),
                                });
                            }}
                            layers={renderer.deckLayers}
                            layerOrder={unifiedLayerOrder}
                            vivLayerProps={
                                renderer.vivLayerProps.length > 0
                                    ? renderer.vivLayerProps
                                    : undefined
                            }
                            deckProps={deckProps}
                        />
                    </div>
                </div>
                {tooltipPortal}
            </SpatialCanvasRendererProvider>
        );
    },
);

function adaptSpatialDataConfig(
    originalConfig: SpatialDataMdvReactConfig,
    dataStore: DataStore,
) {
    const config = { ...scatterDefaults, ...originalConfig };
    if (!dataStore.regions) {
        throw new Error("unexpected attempt to load spatial chart with no regions in datasource");
    }
    config.param = [...dataStore.regions.position_fields];
    if (typeof config.contourParameter === "string") {
        const column = dataStore.columnIndex[config.contourParameter];
        if (!column || allNumeric([column])) {
            config.contourParameter = dataStore.regions.region_field;
        }
    }
    config.viv = applyDefaultChannelState(config.viv);
    if (!config.spatialLayerStack) {
        config.spatialLayerStack = {
            stackOrder: [],
            entries: {},
            spatialLayers: {},
        };
    }
    return config;
}

class SpatialDataMdvReact extends BaseReactChart<SpatialDataMdvReactConfig> {
    declare dataStore: DataStore;
    vivStores: VivContextType;
    spatialCanvasStore: SpatialCanvasStoreApi;
    spatialRendererContext: SpatialCanvasRendererValue | null = null;
    layerDialog?: SpatialLayerDialogReactWrapper;
    ignoreStateUpdate = false;

    get viewerStore() {
        return this.vivStores?.viewerStore;
    }

    constructor(
        dataStore: DataStore,
        div: HTMLDivElement,
        originalConfig: SpatialDataMdvReactConfig,
    ) {
        const config = adaptSpatialDataConfig(originalConfig, dataStore);
        super(dataStore, div, config, SpatialDataChartRoot);
        this.colorByColumn(config.color_by);
        makeObservable(this, {
            colorBy: observable,
            colorByColumn: action,
            colorByDefault: action,
            spatialRendererContext: observable.ref,
        });
        this.vivStores = createVivStores();
        this.spatialCanvasStore = createSpatialCanvasStore();
        if (config.spatialLayerStack) {
            this.syncSpatialLayerStack(config.spatialLayerStack);
        }
        this.addMenuIcon("fas fa-layer-group", "Manage Layers").addEventListener(
            "click",
            () => {
                if (!this.layerDialog) {
                    this.layerDialog = new SpatialLayerDialogReactWrapper(this);
                    this.dialogs.push(this.layerDialog);
                }
            },
        );
    }

    syncSpatialLayerStack(stack: SpatialLayerStackConfig) {
        const { layers, layerOrder } = stackToCanvasState(stack);
        const store = this.spatialCanvasStore.getState();
        store.reset();
        for (const layerId of layerOrder) {
            const layer = layers[layerId];
            if (layer) {
                store.addLayer(layer);
            }
        }
        store.reorderLayers(layerOrder);
    }

    colorBy?: (i: number) => [r: number, g: number, b: number];

    colorByColumn(col?: VivMdvReactConfig["color_by"]) {
        if (!col) return this.colorByDefault();
        this.config.color_by = col;
        //@ts-expect-error legacy color_by options are normalised at runtime by BaseChart.
        this.colorBy = this.getColorFunction(col, true);
    }

    colorByDefault() {
        this.config.color_by = undefined;
        this.colorBy = undefined;
    }

    getColorOptions() {
        return {
            colorby: "all",
        };
    }

    getSettings() {
        const config = this.config;
        const settings = super.getSettings();
        const filters = config.category_filters.map((filter) => {
            const values = this.dataStore.columnIndex[filter.column]?.values?.slice();
            if (!values) throw `failed assertion that we should have a categorical '${filter.column}' here`;
            values.unshift("all");
            return g({
                type: "multidropdown",
                label: `'${filter.column}' filter`,
                current_value: toArray(filter.category),
                values: [values],
                func: (value) => {
                    filter.category = value;
                    config.category_filters = config.category_filters.slice();
                },
            });
        });
        const dataStore = this.dataStore;
        const imageRegionKeys = Object.keys(dataStore.regions?.all_regions ?? {}).filter(
            (regionKey) => dataStore.regions?.all_regions[regionKey].spatial,
        );
        const images = imageRegionKeys.map((regionKey) => ({
            name: regionKey,
            value: regionKey,
        }));

        return settings.concat([
            g({
                type: "dropdown",
                label: `SpatialData (${dataStore.getColumnName(dataStore.regions?.region_field)})`,
                current_value: config.region,
                values: [images, "name", "value"],
                func(value) {
                    if (config.title === config.region) {
                        config.title = value;
                    }
                    config.region = value;
                    config.background_filter.category = value;
                },
            }),
            ...getSharedScatterSettings(config, {
                chart: this,
                includeDensitySettings: true,
                includePointShape: true,
            }),
            g({
                type: "folder",
                label: "Category Filters",
                current_value: filters,
            }),
        ]);
    }

    getConfig() {
        const config = super.getConfig();
        if (this.vivStores) {
            const viewer = this.vivStores.viewerStore.getState();
            config.viv = {
                ...config.viv,
                viewerStore: {
                    viewState: viewer.viewState
                        ? {
                              target: viewer.viewState.target,
                              zoom: viewer.viewState.zoom,
                          }
                        : null,
                },
            };
        }
        if (this.config.spatialLayerStack) {
            config.spatialLayerStack = this.config.spatialLayerStack;
        }
        return config;
    }
}

BaseChart.types.SpatialDataMdvRegionReact = {
    ...BaseChart.types.VivMdvRegionReact,
    init: (config, dataStore, extraConfig) => {
        BaseChart.types.VivMdvRegionReact.init?.(config, dataStore, extraConfig);
        config.type = "SpatialDataMdvRegionReact";
    },
    class: SpatialDataMdvReact,
    name: "SpatialData.js Image Viewer (experimental)",
};

export { SpatialDataMdvReact };
export type SpatialDataMdvReactType = typeof SpatialDataMdvReact;
export default SpatialDataMdvReact;
