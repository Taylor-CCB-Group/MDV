import { useEffect } from "react";
import { useViewerStoreApi } from "./components/avivatorish/state";
import { useChartID } from "./hooks";
import type { VivMDVReact } from "./components/VivMDVReact";

export const useViewStateLink = () => {
    const viewerStore = useViewerStoreApi();
    const id = useChartID();
    
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
        const unsubscribe = viewerStore.subscribe(({ viewState }) => {
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
}
