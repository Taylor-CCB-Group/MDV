import { action, observable, toJS } from "mobx";
import type ChartManager from "./ChartManager";
import type { DataSource } from "./charts";
import _ from "lodash";
import { toPng } from "html-to-image";

/** Per–data-source view config (panel layout and optional highlight state). */
export type ViewDataDataSource = {
    panelWidth?: number;
    layout?: string;
    /** Persisted row indices to highlight when the view is loaded. */
    highlight?: number[];
};

/** View state: dataSources keyed by data source name, initialCharts, optional links. */
export type View = {
    name?: string;
    dataSources: Record<string, ViewDataDataSource>;
    initialCharts: Record<string, unknown[]>;
    links?: unknown[];
    viewImage?: string;
    [key: string]: unknown;
};

export type UpdatedColumns = any;
export type MetaData = any;
export type ChartError = any;

export type State =
    | {
          view: View;
          currentView: string;
          all_views: string[] | null;
          updatedColumns: UpdatedColumns;
          metadata: MetaData;
          chartErrors?: ChartError[];
      }
    | undefined;

const removeImageProp = (state: State) => {
    const cloneState = { ...state };
    if (cloneState?.view?.viewImage) {
        // biome-ignore lint/performance/noDelete: setting to undefined leads to an obscure issue with hasUnsavedChanges
        delete cloneState.view.viewImage;
    }
    return cloneState;
};

class ViewManager {
    @observable accessor current_view = "";
    @observable accessor all_views: string[];
    @observable accessor lastSavedState: State;
    private cm: ChartManager;

    constructor(current_view = "", all_views: string[] = []) {
        this.cm = window.mdv.chartManager;
        this.current_view = current_view;
        this.all_views = all_views;
        // todo: uncomment when we fix state issues
        // setTimeout(() => {
        //     this.lastSavedState = this.cm.getState();
        // }, 1000);
    }

    getCleanPrevState() {
        return removeImageProp(toJS(this.lastSavedState));
    }

    getCleanCurrState() {
        const currState = this.cm.getState();
        return removeImageProp(toJS(currState));
    }

    // Setters
    @action
    setView(view: string) {
        this.current_view = view;
    }

    @action
    setAllViews(all_views: string[]) {
        if (_.isEqual(this.all_views, all_views)) {
            return;
        }
        this.all_views = all_views;
    }

    @action
    setLastSavedState(state: State) {
        this.lastSavedState = state;
    }

    // Check for any unsaved changes and show add view dialog
    @action
    checkAndAddView() {
        // todo: uncomment if else when we fix state issues
        // if (this.hasUnsavedChanges()) {
        //     this.cm.showSaveViewDialog(() => this.cm.showAddViewDialog());
        // } else {
        this.cm.showAddViewDialog();
        // }
    }

    // Check for any unsaved changes and change view
    @action
    checkAndChangeView(view: string, isDelete = false) {
        // todo: uncomment if else when we fix state issues
        // if (this.hasUnsavedChanges() && !isDelete) {
        //     this.cm.showSaveViewDialog(() => {
        //         this.changeView(view);
        //     });
        // } else {
        this.changeView(view);
        // }
    }

    // Change the current view
    @action
    changeView(view: string) {
        // coderabbit had a comment about potential undefined state if there is an error here
        // I think it is somewhat ok for now, but more robust rollback logic is worth considering
        try {
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
                // todo: uncomment when we fix state issues
                // check for anything that might escape our more formal logic
                // setTimeout(() => {
                //     if (this.hasUnsavedChanges()) {
                //         // this happens when changing to a view with gridstack layout
                //         console.warn("Unexpected unsaved changes shortly after changing view");
                //         this.hasUnsavedChanges(true);
                //     } else {
                //         console.log("✅ View changed without unexpected unsaved changes being detected");
                //     }
                // }, 500);
            });
        } catch (error) {
            console.error("error while changing view", error);
        }
    }

    async createImageofView() {
        try {
            // aspect ratio doesn't work properly when the window is resized, commenting for now
            const bounds = this.cm.contentDiv.getBoundingClientRect();
            const aspect = bounds.width / bounds.height;
            const dataUrl = await toPng(this.cm.contentDiv, {
                canvasWidth: 250,
                canvasHeight: 250 / aspect,
            });
            return dataUrl;
        } catch (error) {
            console.error("error while creating image", error);
            return "";
        }
    }

    async getViewDetails() {
        try {
            const viewList: { name: string; image: string }[] = [];
            for (const view of this.all_views) {
                const data = await this.cm.viewLoader(view);
                viewList.push({ name: view, image: data?.viewImage });
            }
            return viewList;
        } catch (error) {
            console.error("error while getting view details", error);
            return [];
        }
    }

    // Save the current state
    @action
    async saveView(errorHandler?: (state: State) => boolean) {
        try {
            const imageUrl = await this.createImageofView();
            const state = this.cm.getState();
            if (state.chartErrors.length > 0) {
                // handling the errors differently if errorHandler is supplied
                if (errorHandler) {
                    const handled = errorHandler(state);
                    // handle it differently if required
                    if (!handled) {
                        throw state.chartErrors;
                    }
                } else {
                    throw state.chartErrors;
                }
            }
            state.view.viewImage = imageUrl;
            this.cm._callListeners("state_saved", state);
            this.setLastSavedState(state);
        } catch (error) {
            console.error("error while saving view", error);
        }
    }

    // Add a new view
    @action
    async addView(viewName: string, checkedDs: { [name: string]: boolean }, isCloneView: boolean) {
        try {
            const { viewData, dsIndex, contentDiv } = this.cm;
            // Add and set the new view only if it doesn't already exist
            if (!this.all_views.includes(viewName)) {
                this.setAllViews([...this.all_views, viewName]);
            }
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
                //add the required datasources (checkedIds) to the view
                for (const name in checkedDs) {
                    if (checkedDs[name]) {
                        state.view.initialCharts[name] = [];
                        state.view.dataSources[name] = {};
                    }
                }

                contentDiv.innerHTML = "";
                await this.cm._init(state.view);
                await this.saveView();
            } else {
                // Get current state before clearing
                const state = this.cm.getState();
                
                // Update the view name in the state
                state.view.name = viewName;
                
                // Clear existing gridstack instances
                for (const ds in viewData.dataSources) {
                    if (viewData.dataSources[ds].layout === "gridstack") {
                        const d = dsIndex[ds];
                        if (!d) continue;
                        this.cm.gridStack.destroy(d);
                    }
                }
                
                // Clear content and reinitialize
                contentDiv.innerHTML = "";
                await this.cm._init(state.view);
                // Save the view
                await this.saveView();
            }
        } catch (error) {
            console.error("error in add view", error);
        }
    }

    // Delete current view
    @action
    deleteView() {
        try {
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
                this.checkAndChangeView(nextView, true);
            } else {
                // no other views exist
                this.cm.removeAllCharts();
                this.cm.viewData = {};
                this.cm.showAddViewDialog();
            }
        } catch (error) {
            console.error("error in delete view", error);
        }
    }

    @action
    async saveAsView(viewName: string) {
        await this.addView(viewName, {}, true);
    }

    checkUnsavedState(action: () => void) {
        if (this.hasUnsavedChanges()) {
            this.cm.showSaveViewDialog(action);
        } else {
            action();
        }
    }

    // Helper to check if the current state differs from the last saved state
    hasUnsavedChanges(verbose = false) {
        if (this.lastSavedState === null) return true;

        const currState = this.cm.getState();

        if (currState.chartErrors && currState.chartErrors.length > 0) {
            console.log("errors while getting state", currState.chartErrors);
            return false;
        }

        const currentState = removeImageProp(toJS(currState));
        const prevState = removeImageProp(toJS(this.lastSavedState));

        if (verbose) {
            console.log("inside currentState", currentState);
            console.log("inside prevState", prevState);
            // doesn't seem to be as useful as hoped for showing what's different...
            // console.log("diff", _.differenceWith([currentState], [prevState], _.isEqual));
            // https://stackoverflow.com/a/48181184/279703 with changes...
            // seems to be managing deep diff of objects subject to further testing
            // we may want more of this in the future... should be a utility function.
            const diff = (obj1: any, obj2: any) => {
                return _.reduce(
                    obj1,
                    (result: any, value: any, key: string) => {
                        if (_.isPlainObject(value)) {
                            const d = diff(value, obj2[key]);
                            if (!_.isEmpty(d)) {
                                result[key] = d;
                            }
                        } else if (!_.isEqual(value, obj2[key])) {
                            // different from SO answer, which wouldn't diff arrays.
                            result[key] = diff(value, obj2[key]);
                        } else {
                            delete result[key];
                        }
                        return result;
                    },
                    {},
                );
            };
            // console.log("diff", diff(currentState, prevState));
        }
        // nb, we should consider updatedColumns, maybe some other things...
        // but not the entire state (e.g. not all_views)
        return !_.isEqual(currentState.view, prevState.view);
    }
}

export default ViewManager;
