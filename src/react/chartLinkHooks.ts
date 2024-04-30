import { useEffect } from "react";
import { useMetadata, useViewerStoreApi } from "./components/avivatorish/state";
import { useChartID } from "./hooks";
import type { VivMDVReact } from "./components/VivMDVReact";

export const useViewStateLink = () => {
    const viewerStore = useViewerStoreApi();
    const id = useChartID();
    const metadata = useMetadata();
    useEffect(() => {
        if (!window.mdv.chartManager?.viewData) return;
        const cm = window.mdv.chartManager;
        const { viewData } = cm; // as of now, this won't change in the lifetime of the component - but hope for interactive link edit soon.
        // 'viewData' is currently just the json metadata for the view, but there could be a UI for manipulating links in it.

        const thisChart = cm.getChart(id) as VivMDVReact;
        // find any "view_state" links in the viewData that include this chart's id in "linked_charts"
        if (!viewData.links) return;
        const vsLinks = viewData.links.filter(l => l.type === "view_state" && l.linked_charts.includes(id));
        if (vsLinks.length === 0) return;
        console.log('found view state link(s)', vsLinks);
        // we want to do something like subscribe to our viewState and push changes to the linked charts
        // but make sure we don't create a circular loop of updates
        const unsubscribe = viewerStore.subscribe(({ viewState }) => {
            thisChart.ignoreStateUpdate = true; //<< add a setting for this, make sure we get the logic right
            const originalZoom = viewState.zoom as number;
            const ourPhysicalSize = metadata.Pixels.PhysicalSizeX;
            const ourUnits = metadata.Pixels.PhysicalSizeXUnit;
            // should we consider viewerStore.useLinkedView?
            // best to be clear about what is and isn't similar to Avivator.
            vsLinks.forEach(link => {
                // as VivMDVReact[] is very much not correct here, we should be checking - and that may not be what we want anyway.
                const otherCharts = link.linked_charts.filter(c => c !== id).map(c => cm.getChart(c)) as VivMDVReact[];
                otherCharts.forEach(c => {
                    if (!c) return; // e.g. if the chart has been removed - ideally the link state would be updated to reflect this.
                    if (c.ignoreStateUpdate) return;
                    // todo - viewState may not be directly compatible with other chart's viewState
                    // so there should be a utility function to convert between them - current attempt is not yet correct, and clutters this code.
                    // might entail some extra garbage collection, making a new object each time. So it goes I guess.
                    const otherMeta = c.vivStores?.viewerStore.getState().metadata;
                    if (!otherMeta) return;
                    const otherPhysicalSize = otherMeta.Pixels.PhysicalSizeX;
                    const otherUnits = otherMeta.Pixels.PhysicalSizeXUnit;
                    if (otherUnits !== ourUnits) throw 'physical size units do not match'; //we could probably convert if this is a common case
                    const zoomRatio = ourPhysicalSize / otherPhysicalSize;
                    const zoom = originalZoom * zoomRatio;
                    // this is not right - target needs to be adjusted so that the same point in the image is centred
                    const target = (viewState.target as [number, number]).map((v, i) => v * zoomRatio);
                    const newViewState = { ...viewState, zoom, target };
                    c.viewerStore?.setState({ viewState: newViewState });
                });
            });
            thisChart.ignoreStateUpdate = false;
        });
        return unsubscribe;
    }, [viewerStore, id]);
}
