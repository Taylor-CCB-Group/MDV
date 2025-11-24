import {
    getDefaultInitialViewState,
    ColorPaletteExtension,
    DetailView,
} from "@vivjs-experimental/viv";
import { observer } from "mobx-react-lite";
import { useMemo, useEffect, useRef } from "react";
import { shallow } from "zustand/shallow";
import { useChartSize, useChartID, useConfig, useRegion } from "../hooks";
import SelectionOverlay from "./SelectionOverlay";
import {
    useLoader,
    type OME_TIFF,
    useViewerStoreApi,
    useChannelsStore,
    useViewerStore,
} from "./avivatorish/state";
import { useViewStateLink } from "../chartLinkHooks";
import { useChart } from "../context";
import { SpatialAnnotationProvider, useSpatialLayers } from "../spatial_context";
import { GeoJsonLayer } from "@deck.gl/layers";
import MDVivViewer, { getVivId } from "./avivatorish/MDVivViewer";
import type { VivRoiConfig } from "./VivMDVReact";
import { useProject } from "@/modules/ProjectContext";
import VivContrastExtension from "@/webgl/VivContrastExtension";
import { trace } from "mobx";
import { useOuterContainer } from "../screen_state";
import type { DeckGLProps, OrbitViewState, OrthographicViewState } from "deck.gl";

export type ViewState = ReturnType<typeof getDefaultInitialViewState>; //<< move this / check if there's an existing type

/** somewhat comparable to avivator `<Viewer />` */
export const VivScatter = () => {
    const chart = useChart();
    return (
        <SpatialAnnotationProvider chart={chart}>
            <Main />
        </SpatialAnnotationProvider>
    );
};

const useJsonLayer = (showJson: boolean) => {
    const id = useChartID();
    const { root } = useProject();
    const { json } = useRegion(); // return type is 'any' and we assume 'json' will be a string - but want that to be different in future.
    const layer_id = `json_${getVivId(`${id}detail-react`)}`;
    const layer = useMemo(() => {
        return json
            ? new GeoJsonLayer({
                  id: layer_id,
                  data: `${root}/${json}`,
                  opacity: 1,
                  filled: true,
                  getFillColor: (f) => [255, 255, 255, 150],
                  getLineColor: (f) => [255, 255, 255, 150],
                  getLineWidth: 2,
                  lineWidthMinPixels: 1,
                  getPointRadius: 10,
                  pickable: true,
                  autoHighlight: true,
                  //@ts-expect-error GeoJson getText: might think about using zod to type/validate this
                  getText: (f) => f.properties.DN,
                  getTextColor: [255, 255, 255, 255],
                  getTextSize: 12,
                  textBackground: true,
                  visible: showJson,
              })
            : null;
    }, [json, showJson, layer_id, root]);
    return layer;
};

const Main = observer(() => {
    // type of this to be sorted - before we accessed ome.data, but maybe this is the 'data'...
    const ome = useLoader() as OME_TIFF["data"]; // useOmeTiff();

    const viewerStore = useViewerStoreApi();
    const [width, height] = useChartSize();
    const id = useChartID();
    const detailId = `${id}detail-react`;
    const outerContainer = useOuterContainer();

    // this isn't updating when we tweak the config...
    const { scatterProps, selectionLayer } = useSpatialLayers();
    const { scatterplotLayer, getTooltip } = scatterProps;
    const { showJson } = useConfig<VivRoiConfig>();
    // passing showJson from here to make use of this being `observer`
    const jsonLayer = useJsonLayer(showJson);

    // maybe more efficient to pick out properties like this... but it's very repetitive/verbose
    const {
        colors,
        contrastLimits,
        channelsVisible,
        selections,
        brightness,
        contrast,
    } = useChannelsStore(
        ({
            colors,
            contrastLimits,
            channelsVisible,
            selections,
            brightness,
            contrast,
        }) => {
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
    const vsRef = useRef<ViewState>();
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
    const extensions = useMemo(
        () => [new ColorPaletteExtension(), new VivContrastExtension()],
        [],
    );
    const detailView = useMemo(
        () =>
            new DetailView({
                id: detailId,
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
    const deckProps: Partial<DeckGLProps> = useMemo(
        () => ({
            getTooltip,
            style: {
                zIndex: "-1",
            },
            //todo figure out why GPU usage is so high (and why commenting and then uncommenting this line fixes it...)
            layers: [jsonLayer, scatterplotLayer, selectionLayer],
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
            }
        }),
        [
            scatterplotLayer,
            selectionLayer,
            jsonLayer,
            id,
            getTooltip,
        ],
    );
    if (!viewState) return <div>Loading...</div>; //this was causing uniforms["sizeScale"] to be NaN, errors in console, no scalebar units...
    if (import.meta.env.DEV) trace();
    return (
        <>
            <SelectionOverlay />
            <MDVivViewer
                outerContainer={outerContainer}
                selectionLayer={selectionLayer}
                views={[detailView]}
                layerProps={[layerConfig]}
                viewStates={[{ ...viewState, id: detailId }]}
                // not really expecting OrbitViewState... yet... but we will for 3D.
                onViewStateChange={(e: {viewState: OrthographicViewState | OrbitViewState}) => {
                    viewerStore.setState({
                        viewState: { ...e.viewState, id: detailId },
                    });
                    if (vsDebugDivRef.current)
                        vsDebugDivRef.current.innerText = JSON.stringify(
                            e.viewState,
                            null,
                            2,
                        );
                }}
                deckProps={deckProps}
            />
        </>
    );
});
