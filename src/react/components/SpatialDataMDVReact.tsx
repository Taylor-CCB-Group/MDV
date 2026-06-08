import { readZarr, type SpatialData } from "@spatialdata/core";
import {
    SpatialCanvasViewer,
    type LayerConfig,
    type ViewState as SpatialCanvasViewState,
} from "@spatialdata/vis";
import type { Layer } from "@deck.gl/core";
import type {
    DeckGLProps,
    OrbitViewState,
    OrthographicViewState,
    PickingInfo,
} from "deck.gl";
import { GeoJsonLayer } from "@deck.gl/layers";
import { observer } from "mobx-react-lite";
import { action, makeObservable, observable } from "mobx";
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
import { useProject } from "@/modules/ProjectContext";
import useGateLayers from "../hooks/useGateLayers";
import { useOuterContainerDeckTooltip } from "../hooks/useOuterContainerDeckTooltip";
import { getSharedScatterSettings } from "./sharedScatterSettings";
import { useOuterContainer } from "../screen_state";
import { scatterDefaults } from "../scatter_state";
import { g, toArray } from "@/lib/utils";

type SpatialRegionMetadata = {
    spatial?: {
        file?: string;
        coordinate_system?: string;
    };
    json?: string;
};

type LoadState =
    | { status: "idle"; spatialData: null; error: null }
    | { status: "loading"; spatialData: null; error: null }
    | { status: "loaded"; spatialData: SpatialData; error: null }
    | { status: "error"; spatialData: null; error: Error };

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

function elementHasCoordinateSystem(element: unknown, coordinateSystem: string) {
    if (!element || typeof element !== "object") return false;
    if (!("coordinateSystems" in element)) return true;
    const systems = element.coordinateSystems;
    return Array.isArray(systems) && systems.includes(coordinateSystem);
}

function firstElementKey(
    elements: Record<string, unknown> | undefined,
    coordinateSystem: string,
) {
    return Object.entries(elements ?? {}).find(([, element]) =>
        elementHasCoordinateSystem(element, coordinateSystem),
    )?.[0] ?? null;
}

function buildSpatialLayerConfig(
    spatialData: SpatialData | null,
    coordinateSystem: string | null,
) {
    if (!spatialData || !coordinateSystem) {
        return { layers: {}, layerOrder: [] };
    }

    const layers: Record<string, LayerConfig> = {};
    const layerOrder: string[] = [];
    const imageKey = firstElementKey(spatialData.images, coordinateSystem);
    if (imageKey) {
        const id = `spatialdata-image-${imageKey}`;
        layers[id] = {
            id,
            type: "image",
            visible: true,
            opacity: 1,
            elementKey: imageKey,
        };
        layerOrder.push(id);
    }

    return { layers, layerOrder };
}

function SpatialDataChartRoot() {
    const { vivStores } = useChart<VivMdvReactConfig, SpatialDataMdvReact>();
    return (
        <VivProvider vivStores={vivStores}>
            <SpatialDataMainChart />
        </VivProvider>
    );
}

const SpatialDataMainChart = observer(() => {
    const chart = useChart<VivMdvReactConfig, SpatialDataMdvReact>();
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

const useJsonLayer = (showJson: boolean) => {
    const id = useChartID();
    const { root } = useProject();
    const region = useRegion();
    const regionJson = getSpatialRegionMetadata(region)?.json;
    const layerId = `json_spatialdata_${id}`;
    return useMemo(() => {
        return regionJson
            ? new GeoJsonLayer({
                  id: layerId,
                  data: `${root}/${regionJson}`,
                  opacity: 1,
                  filled: true,
                  getFillColor: [255, 255, 255, 150],
                  getLineColor: [255, 255, 255, 150],
                  getLineWidth: 2,
                  lineWidthMinPixels: 1,
                  getPointRadius: 10,
                  pickable: true,
                  autoHighlight: true,
                  getText: (feature: { properties?: { DN?: string } }) => feature.properties?.DN,
                  getTextColor: [255, 255, 255, 255],
                  getTextSize: 12,
                  textBackground: true,
                  visible: showJson,
              })
            : null;
    }, [regionJson, showJson, layerId, root]);
};

const SpatialDataViewer = observer(
    ({
        setHoveredField,
    }: {
        setHoveredField: (fieldId: FieldName | null) => void;
    }) => {
        const rawRegion = useRegion();
        const region = getSpatialRegionMetadata(rawRegion);
        const coordinateSystem = region?.spatial?.coordinate_system ?? null;
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
        const { showJson } = useConfig<VivRoiConfig>();
        const jsonLayer = useJsonLayer(showJson);
        const {
            gateLabelLayer,
            gateDisplayLayer,
            controllerOptions,
        } = useGateLayers();
        const config = useConfig<DualContourLegacyConfig>();
        const legendFields = useFieldContourLegend(config.densityFields);
        const showLegend = config.field_legend.display;
        const legendPosition = { x: 10, y: 10 };

        useViewStateLink();

        useEffect(() => {
            if (scatterProps.viewState) {
                viewerStore.setState({ viewState: scatterProps.viewState });
            }
        }, [scatterProps.viewState, viewerStore.setState]);

        const { layers, layerOrder } = useMemo(
            () => buildSpatialLayerConfig(spatialData, coordinateSystem),
            [spatialData, coordinateSystem],
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

        const deckLayers = useMemo(
            () =>
                [
                    jsonLayer,
                    greyScatterplotLayer,
                    scatterplotLayer,
                    gateDisplayLayer,
                    selectionLayer,
                    gateLabelLayer,
                ].filter((layer) => layer !== null) as Layer[],
            [
                gateLabelLayer,
                gateDisplayLayer,
                scatterplotLayer,
                greyScatterplotLayer,
                selectionLayer,
                jsonLayer,
            ],
        );

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
                    <SpatialCanvasViewer
                        spatialData={spatialData}
                        coordinateSystem={coordinateSystem}
                        layers={layers}
                        layerOrder={layerOrder}
                        viewState={spatialViewState}
                        onViewStateChange={(nextViewState) => {
                            viewerStore.setState({
                                viewState: toMdvViewState(
                                    nextViewState,
                                    viewerStore.getState().viewState,
                                ),
                            });
                        }}
                        deckLayers={deckLayers}
                        deckProps={deckProps}
                        renderTooltip={false}
                        showLoadingOverlay
                        autoFit
                        tooltipContainer={outerContainer}
                        style={{
                            width,
                            height,
                        }}
                    />
                </div>
                {tooltipPortal}
            </>
        );
    },
);

function adaptSpatialDataConfig(
    originalConfig: VivMdvReactConfig,
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
    return config;
}

class SpatialDataMdvReact extends BaseReactChart<VivMdvReactConfig> {
    declare dataStore: DataStore;
    vivStores: VivContextType;
    ignoreStateUpdate = false;

    get viewerStore() {
        return this.vivStores?.viewerStore;
    }

    constructor(
        dataStore: DataStore,
        div: HTMLDivElement,
        originalConfig: VivMdvReactConfig,
    ) {
        const config = adaptSpatialDataConfig(originalConfig, dataStore);
        super(dataStore, div, config, SpatialDataChartRoot);
        this.colorByColumn(config.color_by);
        makeObservable(this, {
            colorBy: observable,
            colorByColumn: action,
            colorByDefault: action,
        });
        this.vivStores = createVivStores();
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
                type: "check",
                label: "show json layer",
                current_value: config.showJson || false,
                func: (showJson) => {
                    config.showJson = showJson;
                },
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

export type SpatialDataMdvReactType = typeof SpatialDataMdvReact;
export default 42;
