import { action, observable } from "mobx";
import type ChartManager from "./ChartManager";
import type { DataSource } from "./charts";
import _ from 'lodash';

export type View = any;
export type UpdatedColumns = any;
export type MetaData = any;

export type State = {
    view: View;
    currentView: string;
    all_views: string[] | null;
    updatedColumns: UpdatedColumns;
    metadata: MetaData;
} | undefined;

class ViewManager {
    @observable accessor current_view = "";
    @observable accessor all_views: string[];
    @observable accessor lastSavedState: State;
    private cm: ChartManager;

    constructor(current_view = "", all_views: string[] = []) {
        this.cm = window.mdv.chartManager;
        this.current_view = current_view;
        this.all_views = all_views;
        setTimeout(() => {
            this.lastSavedState = this.cm.getState();
        }, 1000);
    }

    @action
    setView(view: string) {
        this.current_view = view;
    }

    @action
    setAllViews(all_views: string[]) {
        this.all_views = all_views;
    }

    @action
    setLastSavedState(state: State) {
        this.lastSavedState = state;
    }

    @action
    _changeView(view: string) {
        const { viewData, dsIndex, contentDiv } = this.cm;
        for (const ds in this.cm.viewData.dataSources) {
            if (viewData.dataSources[ds].layout === "gridstack") {
                this.cm.gridStack.destroy(dsIndex[ds] as DataSource);
            }
        }
        this.cm.removeAllCharts();
        contentDiv.innerHTML = "";
        this.setView(view);
        this.cm.viewLoader(view).then(async (data: object) => {
            await this.cm._init(data);
            const state = this.cm.getState();
            this.setLastSavedState(state);
        });
    }

    @action
    changeView(view: string, isDelete=false) {
        if (this.hasUnsavedChanges() && !isDelete) {
            this.cm.showSaveViewDialog(() => {
                this._changeView(view);
            });
        } else {
            this._changeView(view);
        }
        
    }

    // Helper to check if the current state differs from the last saved state
    hasUnsavedChanges() {
        const currentState = JSON.parse(JSON.stringify(this.cm.getState()));
        const prevState = JSON.parse(JSON.stringify(this.lastSavedState));
        if (this.lastSavedState === null) return true;
        return !_.isEqual(currentState, prevState);
    }


    @action
    saveView() {
        const state = this.cm.getState();
        this.cm._callListeners("state_saved", state);
        this.setLastSavedState(state);
    }

    @action
    addView(viewName: string, checkedDs: {[name: string]: boolean;}, isCloneView: boolean) {

        const {viewData, dsIndex, contentDiv} = this.cm;
        this.setAllViews([...this.all_views, viewName]);

        // Optionally make it the current view
        this.setView(viewName);
        if (!isCloneView) {
            //remove all charts and links
            for (const ds in viewData.dataSources) {
                if (viewData.dataSources[ds].layout === "gridstack") {
                    const d = dsIndex[ds];
                    if (!d) continue;
                    this.cm.gridStack.destroy(d);
                }
            }
            this.cm.removeAllCharts();
            viewData.links = [];
            const state = this.cm.getState();
            state.view.initialCharts = {};
            state.view.dataSources = {};
            //only one datasource
            if (Object.keys(viewData.dataSources)?.length === 1) {
                const name = Object.keys(viewData.dataSources)?.[0];
                state.view.initialCharts[name] = [];
                state.view.dataSources[name] = {};
            } else {
                for (const ds in dsIndex) {
                    if (checkedDs[ds]) {
                        state.view.initialCharts[ds] = [];
                        state.view.dataSources[ds] = {};
                    }
                }
            }
            contentDiv.innerHTML = "";
            this.cm._init(state.view);
        } else {
            const state = this.cm.getState();
            console.log("state add new: ", state);
            this.cm._init(state.view);
        }
    }

    @action
    deleteView() {
        //remove the view choice and change view to the next one
        const view = this.current_view;

        // update the views
        const updatedViews = this.all_views.filter((v) => v !== view);
        this.setAllViews(updatedViews);

        const state = this.cm.getState();

        //want to delete view and update any listeners
        state.view = null;

        this.cm._callListeners("state_saved", state);

        if (updatedViews.length > 0) {
            // set current view to initial view
            const nextView = updatedViews[0];
            this.setView(nextView);
            this.changeView(nextView);
        } else {
            // no other views exist
            this.cm.removeAllCharts();
            this.cm.viewData = {};
            this.cm.showAddViewDialog();
        }
    }

};

export default ViewManager;