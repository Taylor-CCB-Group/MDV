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
import { getLocalCompressedBinaryDataLoader } from "../dataloaders/DataLoaders.js";

document.addEventListener("DOMContentLoaded", () => loadData());

// if URLSearchParams has a 'dir' parameter, use that as the data directory.
const urlParams = new URLSearchParams(window.location.search);
const dir = urlParams.get('dir') || prompt("Enter data URL", "https://mdvstatic.netlify.app/ytrap2");
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
    let resp = await fetch(url);//, { mode: "no-cors" });
    const config = await resp.json();
    rewriteBaseUrlRecursive(config);
    return config;
}

async function loadData() {
    const datasources = await fetchAndPatchJSON(`${root}/datasources.json`);
    const config = await fetchAndPatchJSON(`${root}/state.json`);
    const views = await fetchAndPatchJSON(`${root}/views.json`);
    const dataLoader = {
        function: getLocalCompressedBinaryDataLoader(datasources, root),
        viewLoader: async (view) => views[view]
    }
    const cm = new ChartManager("app1", datasources, dataLoader, config);

};