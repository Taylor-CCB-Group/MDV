/// <reference types="vite/client" />
import "./all_css";
import ChartManager from "../charts/ChartManager.js";
// import "../charts/RowSummaryBox.js"; //> should this be in ChartManager along with default charts? how useful is it?

import TagModel from '../table/TagModel';
import AnnotationDialog from "../charts/dialogs/AnnotationDialog";
import { BaseDialog } from "../utilities/Dialog";
import SideEffect from "../charts/dialogs/AnnotationDialogReact";
import { getDataLoader } from "../dataloaders/DataLoaderUtil";
console.log(SideEffect);

const flaskURL = "http://localhost:5050";

document.addEventListener("DOMContentLoaded", () => loadData());

// if URLSearchParams has a 'dir' parameter, use that as the data directory.
const urlParams = new URLSearchParams(window.location.search);
// if we're in a popout window, ignore the dir parameter and don't load data
const isPopout = urlParams.get('popout') === "true";
// if there is no dir parameter, use the flaskURL to proxy requests to the python server
const dir = urlParams.get('dir') || (isPopout ? '' : flaskURL);
const root = dir.endsWith("/") ? dir.substring(0, dir.length-1) : dir;
//hack to load data from local API
//TODO - make the API allow serving multiple projects?
const staticFolder = dir !== flaskURL;

// set title of page to the data directory
document.title = staticFolder ? 'MDV - local project' : `MDV - ${root}`;

function rewriteBaseUrlRecursive(config) {
    if (Array.isArray(config)) {
        for (const item of config) {
            rewriteBaseUrlRecursive(item);
        }
        return;
    }
    for (const key in config) {
        if (key === "base_url") {
            config[key] = config[key].replace("./", `${root}/`);
        } else if (typeof config[key] === "object") {
            rewriteBaseUrlRecursive(config[key]);
        }
    }
}

async function fetchAndPatchJSON(url) {
    let resp = await fetch(url)//, { mode: "no-cors" });
    const config = await resp.json();
    rewriteBaseUrlRecursive(config);
    return config;
}

// TODO make a better type for this, put it somewhere more sensible.
export type Datasource = { 
    name: string, 
    columns: { name: string, type: string }[], images?: any, size: number, columnGroups?: any[] 
};


async function loadData() {
    // setupDebug();
    if (isPopout) return;
    const datasources = await fetchAndPatchJSON(`${root}/datasources.json`) as Datasource[];
    const config = await fetchAndPatchJSON(`${root}/state.json`);
    config.popouturl = undefined;
    const views = await fetchAndPatchJSON(`${root}/views.json`);

    const dataLoader = getDataLoader(staticFolder, datasources, views, dir);
    const cm = new ChartManager("app1", datasources, dataLoader, config);
    function extraFeatures(i: number) {
        const dsName = datasources[i].name;
        const ds = cm.dsIndex[dsName];
        const tagModel = new TagModel(ds.dataStore);
        // cm.dsIndex[dsName].menuBar is undefined... so I'm deferring this call.
        // should it be in the viewLoader callback? no ref to cm passed there.
        setTimeout(() => {
            // TODO - add a 'save' button, if supported (another URLSearchParam? Security?)
            cm.addMenuIcon(dsName, "fas fa-tags", "Tag Annotation", () => { new AnnotationDialog(ds.dataStore, tagModel); });
            if (import.meta.env.DEV) {
                cm.addMenuIcon(dsName, "fas fa-tags", "Tag Annotation (react)", () => { new BaseDialog.experiment['AnnotationDialogReact'](ds.dataStore, tagModel); });
            }
            cm.addMenuIcon(dsName, "fas fa-spinner", "Pre-Load Data", async () => { 
                const columns = datasources[i].columns.map(c => c.name);
                cm.loadColumnSet(columns, dsName, () => { console.log("done loadColumnSet"); });
            });
        }, 0);
    }
    datasources.forEach((ds, i) => extraFeatures(i));
};