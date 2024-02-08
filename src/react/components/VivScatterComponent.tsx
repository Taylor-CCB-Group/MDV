import { getDefaultInitialViewState, ColorPaletteExtension, DetailView, VivViewer } from "@hms-dbmi/viv";
import { observer } from "mobx-react-lite";
import { useState, useLayoutEffect, useMemo, useEffect, useRef } from "react";
import { shallow } from "zustand/shallow";
import { useChartSize, useChartID } from "../hooks";
import { useScatterplotLayer } from "../scatter_state";
import SelectionOverlay from "./SelectionOverlay";
import { useLoader, OME_TIFF, useViewerStoreApi, useChannelsStore, useViewerStore } from "./avivatorish/state";
import { VivMDVReact } from "./VivMDVReact";

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
    // I want to make links between viewStates of different charts... this could be moved to a hook.
    useEffect(() => {
        if (!window.mdv.chartManager?.viewData) return;
        const cm = window.mdv.chartManager;
        const { viewData } = cm; // as of now, this won't change in the lifetime of the component - but hope for interactive link edit soon.
        // 'viewData' is currently just the json metadata for the view, but there could be a UI for manipulating links in it.
        
        const thisChart = cm.getChart(id) as VivMDVReact;
        // find any "view_state" links in the viewData that include this chart's id in "linked_charts"
        const vsLinks = viewData.links.filter(l => l.type === "view_state" && l.linked_charts.includes(id));
        if (vsLinks.length === 0) return;
        console.log('found view state link(s)', vsLinks);
        // we want to do something like subscribe to our viewState and push changes to the linked charts
        // but make sure we don't create a circular loop of updates
        thisChart.viewerStore = viewerStore;
        const unsubscribe = viewerStore.subscribe(({viewState}) => {
            thisChart.ignoreStateUpdate = true;
            vsLinks.forEach(link => {
                const otherCharts = link.linked_charts.filter(c => c !== id).map(c => cm.getChart(c));
                otherCharts.forEach(c => {
                    if (c.ignoreStateUpdate) return;
                    c.viewerStore?.setState({ viewState });
                });
            });
            thisChart.ignoreStateUpdate = false;
        });
        return unsubscribe;
    }, [viewerStore, id]);
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
        // _animate: true,
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