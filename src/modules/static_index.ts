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
import { fetchAndPatchJSON, getDataLoader, getPostData, setProjectRoot } from "../dataloaders/DataLoaderUtil";

const flaskURL = "";

document.addEventListener("DOMContentLoaded", () => loadData());

// if URLSearchParams has a 'dir' parameter, use that as the data directory.
const urlParams = new URLSearchParams(window.location.search);
// if we're in a popout window, ignore the dir parameter and don't load data
const isPopout = urlParams.get('popout') === "true";
// if there is no dir parameter, use the flaskURL to proxy requests to the python server
const dir = urlParams.get('dir') || (isPopout ? '' : flaskURL);
const root = dir.endsWith("/") ? dir.substring(0, dir.length-1) : dir;
//hack to load data from local API
const staticFolder = !dir.startsWith("/project");

// set title of page to the data directory
document.title = !staticFolder ? 'MDV - local project' : `MDV - ${root}`;
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
    const datasources = await fetchAndPatchJSON(`${root}/datasources.json`, root) as Datasource[];
    const config = await fetchAndPatchJSON(`${root}/state.json`, root);
    config.popouturl = undefined;
    const views = await fetchAndPatchJSON(`${root}/views.json`, root);

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
    }
    //todo fix searchParams when changing view.
    const cm = new ChartManager("app1", datasources, dataLoader, config, listener);

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