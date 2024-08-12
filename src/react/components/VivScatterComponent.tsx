import { getDefaultInitialViewState, ColorPaletteExtension, DetailView, VivViewer } from "@hms-dbmi/viv";
import { observer } from "mobx-react-lite";
import { useState, useLayoutEffect, useMemo, useEffect, useRef } from "react";
import { shallow } from "zustand/shallow";
import { useChartSize, useChartID, useConfig, useRegion } from "../hooks";
import { useScatterplotLayer } from "../scatter_state";
import SelectionOverlay from "./SelectionOverlay";
import { useLoader, type OME_TIFF, useViewerStoreApi, useChannelsStore, useViewerStore } from "./avivatorish/state";
import { useViewStateLink } from "../chartLinkHooks";
import { useChart } from "../context";
import { SpatialAnnotationProvider, useRange } from "../spatial_context";
import { GeoJsonLayer, PolygonLayer } from "deck.gl/typed";
import MDVivViewer, { getVivId } from "./avivatorish/MDVivViewer";
import type { VivRoiConfig } from "./VivMDVReact";
import { useProject } from "@/modules/ProjectContext";
import VivContrastExtension from "@/webgl/VivContrastExtension";
import { trace } from "mobx";
import { useOuterContainer } from "../screen_state";

export type ViewState = ReturnType<typeof getDefaultInitialViewState>; //<< move this / check if there's an existing type

/** somewhat comparable to avivator `<Viewer />` */
export const VivScatter = () => {
    const chart = useChart();
    return <SpatialAnnotationProvider chart={chart}><Main /></SpatialAnnotationProvider>
}

const useRectLayer = () => {
    const id = useChartID();
    const { start, end } = useRange();
    // note: viv is very picky about layer ids
    const layer_id = `rect_${getVivId(`${id}detail-react`)}`;
    const polygonLayer = useMemo(() => {
        const data = [
            [start, [end[0], start[1]], end, [start[0], end[1]]]
        ];
        return new PolygonLayer({
            id: layer_id,
            data,

            getPolygon: d => d,
            getFillColor: [140, 140, 140, 50],
            getLineColor: [255, 255, 255, 200],
            getLineWidth: 1,
            lineWidthMinPixels: 1,
            // fillOpacity: 0.1, //not working? why is there a prop for it if it doesn't work?
            // opacity: 0.2,
        });
    }, [start, end, layer_id]);
    return polygonLayer;
}

const useJsonLayer = () => {
    const id = useChartID();
    const { root } = useProject();
    const { showJson } = useConfig<VivRoiConfig>();
    const { json } = useRegion(); // return type is 'any' and we assume 'json' will be a string - but want that to be different in future.
    const layer_id = `json_${getVivId(`${id}detail-react`)}`;
    const layer = useMemo(() => {
        return json ? new GeoJsonLayer({
            id: layer_id,
            data: `${root}/${json}`,
            opacity: 0.25,
            filled: true,
            getFillColor: f => [255, 255, 255, 150],
            getLineColor: f => [f.properties.DN, 255, 255, 150],
            getLineWidth: 2,
            lineWidthMinPixels: 1,
            pickable: true,
            autoHighlight: true,
            getText: f => f.properties.DN,
            getTextColor: [255, 255, 255, 255],
            getTextSize: 12,
            textBackground: true,
            visible: showJson,
        }) : null;
    }, [json, showJson, layer_id, root]);
    return layer;
}


const Main = observer(() => {
    // type of this to be sorted - before we accessed ome.data, but maybe this is the 'data'...
    const ome = useLoader() as OME_TIFF['data'];// useOmeTiff();

    const viewerStore = useViewerStoreApi();
    const [width, height] = useChartSize();
    const id = useChartID();
    const detailId = `${id}detail-react`;
    const outerContainer = useOuterContainer();

    //useSpatialLayers()
    const rectLayer = useRectLayer();
    const scatterProps = useScatterplotLayer();
    const {scatterplotLayer, getTooltip} = scatterProps;
    const jsonLayer = useJsonLayer();

    // maybe more efficient to pick out properties like this... but it's very repetitive/verbose
    const {colors, contrastLimits, channelsVisible, selections, brightness, contrast} = useChannelsStore(
        ({ colors, contrastLimits, channelsVisible, selections, brightness, contrast }) => {
            return { colors, contrastLimits, channelsVisible, selections, brightness, contrast }
        },
        shallow
    );

    const viewState = useViewerStore(store => store.viewState);
    useViewStateLink();
    const vsRef = useRef<ViewState>();
    const vsDebugDivRef = useRef<HTMLPreElement>(null);
    
    useEffect(() => {
        if (!ome) return;
        if (!viewState) {
            //WIP <-- c.f. Avivator's useViewerStore() hook
            // setViewState(getDefaultInitialViewState(ome, { width, height }));
            viewerStore.setState({ viewState: getDefaultInitialViewState(ome, { width, height }) });
        }
    }, [ome, width, height, viewState, viewerStore.setState]);
    const extensions = useMemo(() => [new ColorPaletteExtension(), new VivContrastExtension()], []);
    const detailView = useMemo(() => new DetailView({
        id: detailId,
        snapScaleBar: true,
        width, height
    }), [detailId, width, height]);
    useEffect(() => {
        if (scatterProps.viewState) {
            viewerStore.setState({ viewState: scatterProps.viewState });
            // setViewState(scatterProps.viewState);
            vsRef.current = scatterProps.viewState;
        }
    }, [scatterProps.viewState, viewerStore.setState]);
    // TODO get viv working in popouts (seems to be some spurious feature-detection, should be fixed with new version of viv)
    const layerConfig = useMemo(() => ({
        loader: ome,
        selections,
        contrastLimits,
        extensions,
        colors,
        channelsVisible,
        brightness,
        contrast
    }), [ome, selections, contrastLimits, extensions, colors, channelsVisible, brightness, contrast]);
    const deckProps = useMemo(() => ({
        getTooltip,
        style: {
            zIndex: '-1',
        },
        //todo multiple layers, figure out why GPU usage is so high (and why commenting and then uncommenting this line fixes it...)
        layers: [
            jsonLayer,
            scatterplotLayer, rectLayer, 
        ],
        id: `${id}deck`,
        onAfterRender: () => {
            scatterProps.onAfterRender();
        },
        glOptions: {
            preserveDrawingBuffer: true,
        }
    }), [scatterplotLayer, rectLayer, jsonLayer, id, getTooltip, scatterProps.onAfterRender]);
    if (!viewState) return <div>Loading...</div>; //this was causing uniforms["sizeScale"] to be NaN, errors in console, no scalebar units...
    if (import.meta.env.DEV) trace();
    return (
        <>
            <SelectionOverlay {...scatterProps} />
            <MDVivViewer
                outerContainer={outerContainer}
                views={[detailView]}
                layerProps={[layerConfig]}
                viewStates={[{ ...viewState, id: detailId }]}
                onViewStateChange={e => {
                    viewerStore.setState({ viewState: { ...e.viewState, id: detailId } });
                    if (vsDebugDivRef.current) vsDebugDivRef.current.innerText = JSON.stringify(e.viewState, null, 2);
                }}
                deckProps={deckProps}
            />
        </>
    );
});