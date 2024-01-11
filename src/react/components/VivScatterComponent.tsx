import { getDefaultInitialViewState, ColorPaletteExtension, DetailView, VivViewer } from "@hms-dbmi/viv";
import { observer } from "mobx-react-lite";
import { useState, useLayoutEffect, useMemo } from "react";
import { shallow } from "zustand/shallow";
import { useChartSize, useChartID } from "../hooks";
import { useScatterplotLayer } from "../scatter_state";
import SelectionOverlay from "./SelectionOverlay";
import { useLoader, OME_TIFF, useViewerStoreApi, useChannelsStore } from "./avivatorish/state";

/** somewhat comparable to avivator `<Viewer />` */
export const VivScatter = observer(() => {
    // type of this to be sorted - before we accessed ome.data, but maybe this is the 'data'...
    const ome = useLoader() as OME_TIFF['data'];// useOmeTiff();

    const viewerStore = useViewerStoreApi();
    const [width, height] = useChartSize();
    const id = useChartID();
    const detailId = id + 'detail-react';

    const [scatterplotLayer, getTooltip] = useScatterplotLayer();

    const [colors, contrastLimits, channelsVisible, selections] = useChannelsStore(
        store => [
            store.colors,
            store.contrastLimits,
            store.channelsVisible,
            store.selections
        ],
        shallow
    );

    const [viewState, setViewState] = useState<ReturnType<typeof getDefaultInitialViewState>>();
    useLayoutEffect(() => {
        if (!ome) return;
        if (!viewState) {
            //WIP
            setViewState(getDefaultInitialViewState(ome, { width, height }));
        }
    }, [ome]);
    const extensions = useMemo(() => [new ColorPaletteExtension()], []);
    const detailView = useMemo(() => new DetailView({
        id: detailId,
        snapScaleBar: true,
        width, height
    }), [id, width, height]);
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
    }
    if (!viewState) return <div>Loading...</div>; //this was causing uniforms["sizeScale"] to be NaN, errors in console, no scalebar units...
    return (
        <>
            <SelectionOverlay scatterplotLayer={scatterplotLayer} />
            <VivViewer
                views={[detailView]}
                layerProps={[layerConfigX]}
                viewStates={[{ ...viewState, id: detailId }]}
                onViewStateChange={e => {
                    viewerStore.setState({ viewState: { ...e.viewState, id: detailId } });
                }}

                deckProps={deckProps}
            />
        </>
    );
});