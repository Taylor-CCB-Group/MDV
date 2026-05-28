import { getDefaultInitialViewState, ColorPaletteExtension, DetailView } from "@hms-dbmi/viv";
import { observer } from "mobx-react-lite";
import { useMemo, useEffect, useRef, useState, useCallback } from "react";
import { shallow } from "zustand/shallow";
import { useChartSize, useChartID, useConfig } from "../hooks";
import SelectionOverlay from "./SelectionOverlay";
import FieldContourLegend from "./FieldContourLegend";
import { useFieldContourLegend } from "../contour_state";
import type { DualContourLegacyConfig } from "../contour_state";
import { isChartArrayGridMode } from "./chartArrayGridUtils";
import { useVivDensityGridMode, VIV_SCATTER_DECK_KEY } from "../hooks/useVivDensityGridMode";
import { getConfiguredDensityFieldCount } from "./densityGridUtils";
import { useClonedDeckLayersForDeck } from "./chartArrayDeckLayerClones";
import type { MonkeyPatchEditableGeoJsonLayer } from "@/lib/deckMonkeypatch";
import { useDeckMountEpoch } from "../hooks/useDeckMountEpoch";
import type { FieldName } from "@/charts/charts";
import { useLoader, type OME_TIFF, useViewerStoreApi, useChannelsStore, useViewerStore } from "./avivatorish/state";
import { useViewStateLink } from "../chartLinkHooks";
import { useChart } from "../context";
import { SpatialAnnotationProvider, useSpatialLayers } from "../spatial_context";
import MDVivViewer, { type VivViewerWrapperProps } from "./avivatorish/MDVivViewer";
import { useJsonLayer } from "./vivJsonLayer";
import type { VivRoiConfig } from "./VivMDVReact";
import VivContrastExtension from "@/webgl/VivContrastExtension";
import { useOuterContainer } from "../screen_state";
import type { DeckGLProps, OrbitViewState, OrthographicViewState, PickingInfo } from "deck.gl";
import useGateLayers from "../hooks/useGateLayers";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import { useOuterContainerDeckTooltip } from "../hooks/useOuterContainerDeckTooltip";

export type ViewState = ReturnType<typeof getDefaultInitialViewState>; //<< move this / check if there's an existing type

/** somewhat comparable to avivator `<Viewer />` */
export const VivScatter = () => {
    const chart = useChart();
    const [hoveredField, setHoveredField] = useState<FieldName | null>(null);

    return (
        <SpatialAnnotationProvider chart={chart} hoveredFieldId={hoveredField} setHoveredFieldId={setHoveredField}>
            <div className="flex h-full min-h-0 w-full flex-col">
                <Main hoveredField={hoveredField} setHoveredField={setHoveredField} />
            </div>
        </SpatialAnnotationProvider>
    );
};

export type { JsonLayerDeckContext } from "./vivJsonLayer";
export { useJsonLayer } from "./vivJsonLayer";

const Main = observer(
    ({
        hoveredField,
        setHoveredField,
    }: {
        hoveredField: FieldName | null;
        setHoveredField: (fieldId: FieldName | null) => void;
    }) => {
        // type of this to be sorted - before we accessed ome.data, but maybe this is the 'data'...
        const ome = useLoader() as OME_TIFF["data"]; // useOmeTiff();

        const viewerStore = useViewerStoreApi();
        const [width, height] = useChartSize();
        const id = useChartID();
        const detailId = `${id}detail-react`;
        const outerContainer = useOuterContainer();
        const deckContainerRef = useRef<HTMLDivElement | null>(null);

        // this isn't updating when we tweak the config...
        const { scatterProps, selectionLayer } = useSpatialLayers();
        const { scatterplotLayer, greyScatterplotLayer, getTooltip, setScatterKeyboardActive } = scatterProps;
        const { showJson } = useConfig<VivRoiConfig>();
        const roiConfig = useConfig<VivRoiConfig>();
        const contourConfig = useConfig<DualContourLegacyConfig>();
        const showDensityGrid = isChartArrayGridMode({
            chartType: roiConfig.type,
            dimension: roiConfig.dimension,
            layoutMode: contourConfig.density_mode,
            cellCount: getConfiguredDensityFieldCount(contourConfig.densityFields),
        });
        // passing showJson from here to make use of this being `observer`
        const jsonLayer = useJsonLayer(showJson && !showDensityGrid, "overlay");

        const {
            gateLabelLayer,
            gateDisplayLayer,
            controllerOptions,
        } = useGateLayers();

        // Get field contour legend data
        const legendFields = useFieldContourLegend(contourConfig.densityFields);

        // Legend visibility - fixed position in bottom-left
        const showLegend = contourConfig.field_legend.display;

        // Fixed bottom-left position: 10px from left, 10px from bottom
        const legendPosition = { x: 10, y: 10 };

        const handleFieldHover = (fieldId: FieldName | null) => {
            setHoveredField(fieldId);
        };

        // maybe more efficient to pick out properties like this... but it's very repetitive/verbose
        const { colors, contrastLimits, channelsVisible, selections, brightness, contrast } = useChannelsStore(
            ({ colors, contrastLimits, channelsVisible, selections, brightness, contrast }) => {
                return {
                    colors,
                    contrastLimits,
                    channelsVisible,
                    selections,
                    brightness,
                    contrast,
                };
            },
            shallow,
        );

        const viewState = useViewerStore((store) => store.viewState);
        useViewStateLink();
        const vsRef = useRef<ViewState | null>(null);
        const vsDebugDivRef = useRef<HTMLPreElement>(null);

        useEffect(() => {
            if (!ome) return;
            if (!viewState) {
                //WIP <-- c.f. Avivator's useViewerStore() hook
                // setViewState(getDefaultInitialViewState(ome, { width, height }));
                viewerStore.setState({
                    viewState: getDefaultInitialViewState(ome, { width, height }),
                });
            }
        }, [ome, width, height, viewState, viewerStore.setState]);
        const extensions = useMemo(() => [new ColorPaletteExtension(), new VivContrastExtension()], []);
        const detailView = useMemo(
            () =>
                new DetailView({
                    id: detailId,
                    // @ts-expect-error viv runtime supports this, types do not
                    snapScaleBar: true,
                    width,
                    height,
                }),
            [detailId, width, height],
        );
        useEffect(() => {
            if (scatterProps.viewState) {
                viewerStore.setState({ viewState: scatterProps.viewState });
                // setViewState(scatterProps.viewState);
                vsRef.current = scatterProps.viewState;
            }
        }, [scatterProps.viewState, viewerStore.setState]);
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
                [
                ome,
                selections,
                contrastLimits,
                extensions,
                colors,
                channelsVisible,
                brightness,
                contrast,
            ],
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

        const deckMountEpoch = useDeckMountEpoch(showDensityGrid);

        const overlayDeckLayers = useClonedDeckLayersForDeck(
            [
                jsonLayer,
                greyScatterplotLayer,
                scatterplotLayer,
                gateDisplayLayer,
                selectionLayer,
                gateLabelLayer,
            ],
            deckMountEpoch,
            { mode: "overlay" },
            !showDensityGrid,
        );

        const overlaySelectionLayer = useMemo(
            () =>
                overlayDeckLayers.find((layer) => layer.id.startsWith("selection_")) as
                    | MonkeyPatchEditableGeoJsonLayer
                    | undefined,
            [overlayDeckLayers],
        );

        const deckProps: Partial<DeckGLProps> = useMemo(
            () => ({
                getTooltip: getPortalTooltip,
                layers: overlayDeckLayers,
                id: `${id}deck`,
                // deviceProps: {
                //     webgl: {                    
                //         depth: true,
                //         preserveDrawingBuffer: true,
                //         antialias: true,
                //     },
                // },
                controller: {
                    doubleClickZoom: false,
                    dragPan: controllerOptions.dragPan,
                },
                // deviceProps: {
                //     // todo - get this working more usefully.
                //     debugSpectorJS: true,
                // }
            }),
            [overlayDeckLayers, id, getPortalTooltip, controllerOptions],
        );
        const gridMode = useVivDensityGridMode(showDensityGrid, deckContainerRef);

        if (!viewState) return <div>Loading...</div>; //this was causing uniforms["sizeScale"] to be NaN, errors in console, no scalebar units...

        if (showDensityGrid && gridMode?.blocking) {
            return (
                <>
                    <SelectionOverlay />
                    <div
                        className="min-h-0 min-w-0 flex-1 flex-col"
                        style={{ width: "100%", height: "100%", outline: "none" }}
                    >
                        {gridMode.blocking}
                    </div>
                </>
            );
        }

        const overlayOnViewStateChange = (e: {
            viewState: OrthographicViewState | OrbitViewState;
        }) => {
            viewerStore.setState({
                viewState: { ...e.viewState, id: detailId },
            });
            if (vsDebugDivRef.current) {
                vsDebugDivRef.current.innerText = JSON.stringify(e.viewState, null, 2);
            }
        };

        const activeHandlers = showDensityGrid && gridMode
            ? gridMode.containerHandlers
            : {
                  onPointerDown: suppressTooltipUntilPointerUp,
                  onMouseDown: () => setScatterKeyboardActive(true),
                  onMouseEnter: () => setScatterKeyboardActive(true),
                  onMouseLeave: () => {
                      clearTooltip();
                      setScatterKeyboardActive(false);
                  },
              };

        const activeTooltipPortal = showDensityGrid && gridMode ? gridMode.tooltipPortal : tooltipPortal;
        const gridViewerProps =
            showDensityGrid && gridMode?.viewer
                ? (gridMode.viewer as unknown as VivViewerWrapperProps)
                : null;

        return (
            <>
                <SelectionOverlay />
                {!showDensityGrid && showLegend && legendFields.length > 0 && (
                    <FieldContourLegend
                        fields={legendFields}
                        position={legendPosition}
                        onFieldHover={handleFieldHover}
                    />
                )}
                <div
                    ref={deckContainerRef}
                    className="relative flex min-h-0 min-w-0 flex-1 flex-col"
                    aria-label={showDensityGrid ? "Spatial scatter density grid" : "Spatial scatter plot"}
                    style={{ width: "100%", height: "100%", outline: "none" }}
                    {...activeHandlers}
                >
                    {gridViewerProps ? (
                        <div className="pointer-events-auto absolute inset-0 z-0 overflow-hidden">
                            <MDVivViewer key={VIV_SCATTER_DECK_KEY} useDevicePixels {...gridViewerProps} />
                        </div>
                    ) : (
                        <div className="pointer-events-auto absolute inset-0 z-0 overflow-hidden">
                            <MDVivViewer
                                key={VIV_SCATTER_DECK_KEY}
                                outerContainer={outerContainer}
                                selectionLayer={overlaySelectionLayer ?? selectionLayer}
                                views={[detailView]}
                                layerProps={[layerConfig]}
                                viewStates={[{ ...viewState, id: detailId }]}
                                useDevicePixels={true}
                                onViewStateChange={overlayOnViewStateChange}
                                deckProps={deckProps}
                            />
                        </div>
                    )}
                    {showDensityGrid && gridMode?.layout ? (
                        <div className="relative z-[1] flex min-h-0 min-w-0 flex-1 flex-col pointer-events-none">
                            {gridMode.layout}
                        </div>
                    ) : null}
                </div>
                {activeTooltipPortal}
            </>
        );
});
