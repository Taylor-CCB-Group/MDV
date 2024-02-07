/// <reference types="vite/client" />
import "./all_css";
import "../charts/VivScatterPlotNew";
import HmrHack from "../react/HmrHack";
HmrHack();
import ChartManager from "../charts/ChartManager.js";
// import "../charts/RowSummaryBox.js"; //> should this be in ChartManager along with default charts? how useful is it?

import TagModel from '../table/TagModel';
import AnnotationDialog from "../charts/dialogs/AnnotationDialog";
import { BaseDialog } from "../utilities/Dialog";
import { fetchJsonConfig, getDataLoader, getPostData, setProjectRoot } from "../dataloaders/DataLoaderUtil";
import { changeURLParam } from "./desktop_index";
import BaseChart from "../charts/BaseChart";

// see also basic_index.js for some global mdv stuff... only interested in chartManager for now.
declare global {
    interface Window {
        mdv: {
            ChartManager: typeof ChartManager,
            chartManager?: ChartManager,
            chartTypes?: any,
        };
    }
}
window.mdv = {
    // does this not want to be a singleton - rather than referring to the class?
    // not sure what I'd use the class for - would generally import in a module
    ChartManager,
    chartTypes: BaseChart.types,
}


// todo: change so that *only* ?dir=... is used to determine the root?
// const flaskURL = window.location.pathname;
const { href } = window.location;
const flaskURL = new URL(href).origin;//href.substring(href.indexOf("/project"));


document.addEventListener("DOMContentLoaded", async () => {
    try {
        await loadData();
    } catch (e) {
        alert(`Error loading data: ${e}`)
    }
});

/// --- this section is a bit of a mess, but it works for now --- 
//(copilot suggesting "this is a bit of a mess" is a bit rude, but in this case it's right)
// if URLSearchParams has a 'dir' parameter, use that as the data directory.
const urlParams = new URLSearchParams(window.location.search);
// if we're in a popout window, ignore the dir parameter and don't load data
const isPopout = urlParams.get('popout') === "true";
// if there is no dir parameter, use the flaskURL to proxy requests to the python server
const dir = urlParams.get('dir') || (isPopout ? '' : flaskURL);
function getRoot(dir: string) {
    // const url = new URL(dir);
    // return url.origin + url.pathname;
    return dir.endsWith("/") ? dir.substring(0, dir.length-1) : dir;
}
const root = getRoot(dir);
//hack to load data from local API... TODO XXX make less hacky, rather than more...
//this was sort-of-working as of this writing for `MDVProject.serve()`, as long as the default port 5050 is used...
//need to revisit and actually come up with a proper design.
const staticFolder = false; //!dir.startsWith("/project") && !(window.location.port === "5050") && !dir.endsWith("5050");
const project_name = dir.split("/").pop();
/// --- end of messy section ---

// set title of page to the data directory
document.title = `MDV - ${project_name}`;
if (isPopout) document.title = "MDV popout";


// TODO make a better type for this, put it somewhere more sensible.
export type Datasource = { 
    name: string, 
    columns: { name: string, type: string }[], images?: any, size: number, columnGroups?: any[] 
};


async function loadData() {
    // setupDebug();
    if (isPopout) return;
    setProjectRoot(root);
    // move more of this into DataLoaderUtil (which might get renamed)?
    const datasources = await fetchJsonConfig(`${root}/datasources.json`, root) as Datasource[];
    const config = await fetchJsonConfig(`${root}/state.json`, root);
    config.popouturl = undefined;
    const views = await fetchJsonConfig(`${root}/views.json`, root);
    //is view in the URL
    const view = urlParams.get("view");
    if (config.all_views && view && config.all_views.indexOf(view) !== -1) {
        config["initial_view"] = view;
    }


    const dataLoader = getDataLoader(staticFolder, datasources, views, dir);

    const listener = async (type: string, cm: ChartManager, data: any) => {
        if (type === "state_saved" && !staticFolder) {
            const resp = await getPostData(root+'/save_state', data);
            if (resp.success) {
                cm.createInfoAlert("State saved", {duration: 2000});
                cm.setAllColumnsClean();
            } else {
                cm.createInfoAlert("State save failed", {duration: 3000, type: "danger"});
            }
        }
        if (type === "view_loaded") {
            changeURLParam("view", cm.currentView)
        }
    }
    // this will assign itself to `window.mdv.chartManager` in constructor
    const cm = new ChartManager("holder", datasources, dataLoader, config, listener);

    function extraFeatures(i: number) {
        const dsName = datasources[i].name;
        const ds = cm.dsIndex[dsName];
        const tagModel = new TagModel(ds.dataStore);
        // cm.dsIndex[dsName].menuBar is undefined... so I'm deferring this call.
        // should it be in the viewLoader callback? no ref to cm passed there.
        setTimeout(() => {
            cm.addMenuIcon(dsName, "fas fa-tags", "Tag Annotation", () => { new AnnotationDialog(ds.dataStore, tagModel); });
            if (import.meta.env.DEV) {
                cm.addMenuIcon(dsName, "fas fa-tags", "Tag Annotation (react)", () => { new BaseDialog.experiment['AnnotationDialogReact'](ds.dataStore, tagModel); });
            }
            cm.addMenuIcon(dsName, "fas fa-spinner", "Pre-Load Data", async () => { 
                const columns = datasources[i].columns.map(c => c.name);
                cm.loadColumnSet(columns, dsName, () => { console.log("done loadColumnSet"); });
            });
        }, 100);
    }
    datasources.forEach((ds, i) => extraFeatures(i));
};