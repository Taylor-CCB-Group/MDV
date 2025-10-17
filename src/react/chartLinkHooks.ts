import { useEffect, useId, useMemo, useState } from "react";
import { useMetadata, useViewerStoreApi } from "./components/avivatorish/state";
import { useChartID, useDataSources } from "./hooks";
import type { VivMDVReact } from "./components/VivMDVReact";
import { useDataStore } from "./context";
import type { DataColumn, DataType } from "@/charts/charts";
import { getRowsAsColumnsLinks } from "@/links/link_utils";

export type ChartID = string;
export type ViewStateLink = {
    type: "view_state";
    linked_charts: ChartID[];
}

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
        const vsLinks = viewData.links.filter(
            (l: ViewStateLink) => l.type === "view_state" && l.linked_charts.includes(id),
        );
        if (vsLinks.length === 0) return;
        console.log("found view state link(s)", vsLinks);
        // we want to do something like subscribe to our viewState and push changes to the linked charts
        // but make sure we don't create a circular loop of updates
        const unsubscribe = viewerStore.subscribe(({ viewState }) => {
            thisChart.ignoreStateUpdate = true; //<< add a setting for this, make sure we get the logic right
            const originalZoom = viewState.zoom as number;
            //@ts-expect-error - TODO: fix avivatorish metadata type
            const ourPhysicalSize = metadata?.Pixels.PhysicalSizeX 
            if (!ourPhysicalSize) {
                //! todo - allow this link to work without physical size/viv context.
                console.warn("no physical size in metadata, unexpected metadata format, or used outside viv context?");
                return;
            }
            //@ts-expect-error - TODO: fix avivatorish metadata type
            const ourUnits = metadata.Pixels.PhysicalSizeXUnit;
            // should we consider viewerStore.useLinkedView?
            // best to be clear about what is and isn't similar to Avivator.
            vsLinks.forEach((link: ViewStateLink) => {
                // as VivMDVReact[] is very much not correct here, we should be checking - and that may not be what we want anyway.
                const otherCharts = link.linked_charts
                    .filter((c) => c !== id)
                    .map((c) => cm.getChart(c)) as VivMDVReact[];
                otherCharts.forEach((c) => {
                    if (!c) return; // e.g. if the chart has been removed - ideally the link state would be updated to reflect this.
                    if (c.ignoreStateUpdate) return;
                    // todo - viewState may not be directly compatible with other chart's viewState
                    // so there should be a utility function to convert between them - current attempt is not yet correct, and clutters this code.
                    // might entail some extra garbage collection, making a new object each time. So it goes I guess.
                    const otherMeta =
                        c.vivStores?.viewerStore.getState().metadata;
                    if (!otherMeta) return;
                    //@ts-expect-error - TODO: fix avivatorish metadata type
                    const otherPhysicalSize = otherMeta.Pixels.PhysicalSizeX;
                    if (!otherPhysicalSize) return;
                    //@ts-expect-error - TODO: fix avivatorish metadata type
                    const otherUnits = otherMeta.Pixels.PhysicalSizeXUnit;
                    if (otherUnits !== ourUnits)
                        throw "physical size units do not match"; //we could probably convert if this is a common case
                    const zoomRatio = ourPhysicalSize / otherPhysicalSize;
                    const zoom = originalZoom * zoomRatio;
                    // this is not right - target needs to be adjusted so that the same point in the image is centred
                    const target = (viewState.target as [number, number]).map(
                        (v, i) => v * zoomRatio,
                    );
                    const newViewState = { ...viewState, zoom, target };
                    c.viewerStore?.setState({ viewState: newViewState });
                });
            });
            thisChart.ignoreStateUpdate = false;
        });
        return unsubscribe;
    }, [
        viewerStore,
        id,
        //@ts-expect-error - TODO: fix avivatorish metadata type
        metadata?.Pixels.PhysicalSizeX,
        //@ts-expect-error - TODO: fix avivatorish metadata type
        metadata?.Pixels.PhysicalSizeXUnit,
    ]);
};

export type RowsAsColslink = {
    name_column: string;
    name: string;
    subgroups: {
        [sgName: string]: {
            name: string;
            label: string;
            type: string;
        } 
    }
}

/** returns information about any `rows_as_columns_link`s that exist in the context `dataSource` */
export function useRowsAsColumnsLinks() {
    //- we should have useDataSource() which would work in dialogs not associated with a particular chart.
    //(test in AddChartDialog)
    const dataStore = useDataStore(); //! this whole dataSource vs dataStore thing still confuses me
    if (!dataStore) {
        throw "no dataStore!!!";
    }
    // if it was possible for user to edit this at runtime, we'd want this to be reactive
    //nb, not making this async because of how it's used for deserisation...
    //! we could check initPromise of each link & only return resolved ones...
    return getRowsAsColumnsLinks(dataStore);
}

/** What if the linkIndex might be able to change at runtime?
 * 
 * We could also consider having distinct hooks for highlighted and filtered rows.
 * @returns the text values and indices of the highlighted/filtered rows in a linked dataSource
 */
export function useHighlightedForeignRows(linkIndex = 0) {
    const id = useId();
    const racLink = useRowsAsColumnsLinks()[linkIndex];
    const values = racLink?.link.observableFields || [];
    return values; //maybe don't useState version of this but return the mobx observable directly
}
/** design of this will need to change to account for n-links
 * 
 * Will return an array of DataColumn objects, with loaded data, updated as the highlighted/filtered rows change
 * in the linked DataSource.
 * 
 * We could also consider 
 * - having distinct hooks for highlighted and filtered rows (potentially ways of composing custom filter graphs).
 * - controls for pagination.
 * - returning column objects that are not loaded yet, but will be loaded when they are needed.
 * 
 * 
 * @param max - maximum number of columns to return
 * @returns an array of DataColumn objects, with loaded data, for virtual columns 
 * corresponding to the highlighted/filtered rows in the linked dataSource.
 */
export function useHighlightedForeignRowsAsColumns(max = 10, filter = "", linkIndex = 0) {
    const cols = useHighlightedForeignRows(linkIndex); //not actual cols
    //would like not to have this here - might have some more logic in above hook
    //in particular want to redesign the fieldName being what determines the column
    const racLink = useRowsAsColumnsLinks()[linkIndex];
    //! this next line may be problematic - useMemo is important
    const { link } = useMemo(() => racLink || { link: null }, [racLink]);
    const [columns, setColumns] = useState<DataColumn<DataType>[]>([]);
    const ds = useDataStore();
    useEffect(() => {
        // todo - check whether checking link for null here means we can clean up ForeignRows component
        // (which currently breaks rules of hooks)
        if (!link) {
            return;
        }
        if (cols.length === 0) {
            setColumns([]);
            return;
        }
        const cm = window.mdv.chartManager;
        // c.value & c.index are from the DataStore listener event, now in cols
        // this needs to change for multiple subgroups
        const sg = Object.keys(link.subgroups)[0];
        const f = filter.toLowerCase();
        // todo: consider pagination... pass in a page number, return information about total number of columns etc.
        // moving ds.addColumnsFromFields() up the stack...
        const c = cols.filter(({name}) => name.toLowerCase().includes(f)).slice(0, max).map(({column}) => column);
        // don't setColumns until we have the data...
        // alternatively, they could be lazy-loaded by consuming components.
        // that could be implemented relatively easily, but would lack some batching of requests,
        // which may or may not be useful.
        const cFields = c.map(col => col.field);
        cm.loadColumnSet(cFields, ds.name, () => {
            setColumns(c);
        });
        return; //could consider cancelling any pending requests...
    }, [cols, max, link, ds, filter]);
    return columns;
}
