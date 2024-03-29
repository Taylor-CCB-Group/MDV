/// <reference types="vite/client" />
import 'nouislider/dist/nouislider.min.css'
import "microtip/microtip.css";
import "../../src/css/fontawesome-5.15.3/all.css";
import "../../src/utilities/css/ContextMenu.css";
import "../../src/charts/css/charts.css";
import "../../src/webgl/css/wgl2di.css";
import "../../src/table/css/slickgrid.css";
import ChartManager from "../charts/ChartManager.js";
import "../charts/RowSummaryBox.js";
import "../charts/ImageTableChart.js";
import "../charts/ImageScatterChart.js";
import "../charts/WordCloudChart.js";
import { getLocalCompressedBinaryDataLoader } from "../dataloaders/DataLoaders.js";

import TagModel from '../table/TagModel';
import AnnotationDialog from "../charts/dialogs/AnnotationDialog";
import { BaseDialog } from "../utilities/Dialog";
import SideEffect from "../charts/dialogs/AnnotationDialogReact";
console.log(SideEffect);


document.addEventListener("DOMContentLoaded", () => loadData());

// if URLSearchParams has a 'dir' parameter, use that as the data directory.
const urlParams = new URLSearchParams(window.location.search);
// if we're in a popout window, ignore the dir parameter and don't load data
const isPopout = urlParams.get('popout') === "true";
const dir = urlParams.get('dir') || (isPopout ? '' : prompt("Enter data URL", "https://mdvstatic.netlify.app/ytrap2"));
const root = dir.endsWith("/") ? dir.substring(0, dir.length-1) : dir;
// set title of page to the data directory
document.title = `MDV - ${root}`;

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

async function loadData() {
    if (isPopout) return;
    const datasources = await fetchAndPatchJSON(`${root}/datasources.json`);
    const config = await fetchAndPatchJSON(`${root}/state.json`);
    const views = await fetchAndPatchJSON(`${root}/views.json`);
    // TODO: add a way of specifying a different type of data loader.
    /// perhaps via another URLSearchParam, to avoid messing with datasources spec.
    const dataLoader = {
        function: getLocalCompressedBinaryDataLoader(datasources, root),
        viewLoader: async (view) => views[view]
    }
    const cm = new ChartManager("app1", datasources, dataLoader, config);
    const dsName = datasources[0].name;
    const ds = cm.dsIndex[dsName];
    const tadModel = new TagModel(ds.dataStore);
    // cm.dsIndex[dsName].menuBar is undefined... so I'm deferring this call.
    setTimeout(() => {
        // TODO - add a 'save' button, if supported (another URLSearchParam? Security?)
        cm.addMenuIcon(dsName, "fas fa-tags", "Tag Annotation", () => { new AnnotationDialog(ds.dataStore, tadModel); });
        if (import.meta.env.DEV) {
            cm.addMenuIcon(dsName, "fas fa-tags", "Tag Annotation (react)", () => { new BaseDialog.experiment['AnnotationDialogReact'](ds.dataStore, tadModel); });
        }
        cm.addMenuIcon(dsName, "fas fa-spinner", "Pre-Load Data", async () => { 
            const columns = datasources[0].columns.map(c => c.name);
            cm.loadColumnSet(columns, dsName, () => { console.log("done loadColumnSet"); });
        });
    }, 0);
};