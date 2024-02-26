import { getDefaultInitialViewState, ColorPaletteExtension, DetailView, VivViewer } from "@hms-dbmi/viv";
import { observer } from "mobx-react-lite";
import { useState, useLayoutEffect, useMemo, useEffect, useRef } from "react";
import { shallow } from "zustand/shallow";
import { useChartSize, useChartID } from "../hooks";
import { useScatterplotLayer } from "../scatter_state";
import SelectionOverlay from "./SelectionOverlay";
import { useLoader, OME_TIFF, useViewerStoreApi, useChannelsStore, useViewerStore } from "./avivatorish/state";
import { useViewStateLink } from "../chartLinkHooks";

export type ViewState = ReturnType<typeof getDefaultInitialViewState>; //<< move this / check if there's an existing type

/** somewhat comparable to avivator `<Viewer />` */
export const VivScatter = observer(() => {
    // type of this to be sorted - before we accessed ome.data, but maybe this is the 'data'...
    const ome = useLoader() as OME_TIFF['data'];// useOmeTiff();

    const viewerStore = useViewerStoreApi();
    const [width, height] = useChartSize();
    const id = useChartID();
    const detailId = id + 'detail-react';

    const scatterProps = useScatterplotLayer();
    const {scatterplotLayer, getTooltip} = scatterProps;

    const [colors, contrastLimits, channelsVisible, selections] = useChannelsStore(
        store => [
            store.colors,
            store.contrastLimits,
            store.channelsVisible,
            store.selections
        ],
        shallow
    );

    const viewState = useViewerStore(store => store.viewState);
    useViewStateLink();
    const vsRef = useRef<ViewState>();
    const vsDebugDivRef = useRef<HTMLPreElement>(null);
    
    useLayoutEffect(() => {
        if (!ome) return;
        if (!viewState) {
            //WIP <-- c.f. Avivator's useViewerStore() hook
            // setViewState(getDefaultInitialViewState(ome, { width, height }));
            viewerStore.setState({ viewState: getDefaultInitialViewState(ome, { width, height }) });
        }
    }, [ome]);
    const extensions = useMemo(() => [new ColorPaletteExtension()], []);
    const detailView = useMemo(() => new DetailView({
        id: detailId,
        snapScaleBar: true,
        width, height
    }), [id, width, height]);
    useEffect(() => {
        if (scatterProps.viewState) {
            viewerStore.setState({ viewState: scatterProps.viewState });
            // setViewState(scatterProps.viewState);
            vsRef.current = scatterProps.viewState;
        }
    }, [scatterProps.viewState])
    // TODO get viv working in popouts (not a react thing - happens elsewhere
    // - probably need to handle lost gl context)
    const layerConfigX = {
        loader: ome,
        selections,
        contrastLimits,
        extensions,
        colors,
        channelsVisible,
    }
    const deckProps = {
        getTooltip,
        style: {
            zIndex: '-1',
        },
        layers: [scatterplotLayer],
        id: id + 'deck',
        onAfterRender: () => {
            scatterProps.onAfterRender();
        }
    }
    if (!viewState) return <div>Loading...</div>; //this was causing uniforms["sizeScale"] to be NaN, errors in console, no scalebar units...
    return (
        <>
            <SelectionOverlay {...scatterProps} />
            <pre ref={vsDebugDivRef} 
            style={{ 
                display: 'none', // <-- set to block to debug viewState
                position: 'absolute', bottom: 2, left: 2, zIndex: 1000, backgroundColor: 'rgba(40,40,40,0.5)',
                backdropFilter: 'blur(2px)',
                outline: '1px solid white',
                fontSize: '9px',
                //  pointerEvents: 'none'
            }}
            onClick={() => {
                // center on image
                const vs = getDefaultInitialViewState(ome, { width, height }) as any;
                vs.zoom = (vs.zoom || -10) + Math.random()*0.0001; //salt the zoom to force a re-render
                vsRef.current = vs;
                // setViewState(vs);
                viewerStore.setState({ viewState: vs });
            }}
            >
            </pre>
            <VivViewer
                views={[detailView]}
                layerProps={[layerConfigX]}
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