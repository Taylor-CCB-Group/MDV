/// <reference types="vite/client" />
import "./all_css";
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
import DebugJsonReactWrapper from "@/react/components/DebugJsonDialogReactWrapper";

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




document.addEventListener("DOMContentLoaded", async () => {
    try {
        await loadData();
    } catch (e) {
        // alert(`Error loading data: ${e}`)
        document.head.title = "MDV - Error loading data"; // not actually appearing, hey-ho.
        document.body.style.textAlign = "center";
        document.body.innerHTML = `<h1>Error loading data</h1><p>${e}</p>`;
    }
});

/// --- this section is a bit of a mess, but it works for now ---
//(copilot suggesting "this is a bit of a mess" is a bit rude, but in this case it's an understatement)
// todo: change so that *only* ?dir=... is used to determine the root?
// const flaskURL = window.location.pathname;
const { origin, pathname } = window.location;
const flaskURL = origin + pathname;//new URL(href).origin;//href.substring(href.indexOf("/project"));
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
// - this means some online projects are currently broken, I should fix that. <--
// as of this writing, they work with a `?static` parameter, but that's not a good solution.
const staticFolder = urlParams.get('static') !== null; //!dir.startsWith("/project") && !(window.location.port === "5050") && !dir.endsWith("5050");
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

    // add a button for debugging datasources & views metadata
    // this could be in ChartManager instead, this is convenient for now so we have `root` available.
    // would be better if it appeared with other entry-points, and also having it here means HMR doesn't work.
    const tiptext = "View datasource metadata";
    const debugButton = cm.addMenuIcon('_main', "fas fa-bug", tiptext, async () => {
        const datasources = await fetchJsonConfig(`${root}/datasources.json`, root);
        const views = await fetchJsonConfig(`${root}/views.json`, root);
        
        new DebugJsonReactWrapper({datasources, views});
    });
    debugButton.style.float = "right";
    debugButton.setAttribute("data-microtip-position", "bottom-left");

    function extraFeatures(i: number) {
        const dsName = datasources[i].name;
        const ds = cm.dsIndex[dsName];
        let tagModel: TagModel; // = new TagModel(ds.dataStore);
        function getTagModel() {
            if (!tagModel) tagModel = new TagModel(ds.dataStore);
            return tagModel;
        }
        // ds.dataStore.removeColumn('__tags');
        // cm.dsIndex[dsName].menuBar is undefined... so I'm deferring this call.
        // should it be in the viewLoader callback? no ref to cm passed there.
        setTimeout(() => {
            cm.addMenuIcon(dsName, "fas fa-tags", "Tag Annotation", () => { new AnnotationDialog(ds.dataStore, getTagModel()); });
            if (import.meta.env.DEV) {
                cm.addMenuIcon(dsName, "fas fa-tags", "Tag Annotation (react)", () => { new BaseDialog.experiment['AnnotationDialogReact'](ds.dataStore, getTagModel()); });
            }
            cm.addMenuIcon(dsName, "fas fa-spinner", "Pre-Load Data", async () => {
                const columns = datasources[i].columns.map(c => c.name);
                cm.loadColumnSet(columns, dsName, () => { console.log("done loadColumnSet"); });
            });
        }, 100);
    }
    datasources.forEach((ds, i) => extraFeatures(i));
};
