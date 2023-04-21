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
const dir = '/static/' + urlParams.get('dir') || '/data/ytrap';//http://localhost:8082/';

function rewriteBaseUrlRecursive(config) {
    if (Array.isArray(config)) {
        for (const item of config) {
            rewriteBaseUrlRecursive(item);
        }
        return;
    }
    for (const key in config) {
        if (key === "base_url") {
            config[key] = config[key].replace("./", `${dir}/`);
        } else if (typeof config[key] === "object") {
            rewriteBaseUrlRecursive(config[key]);
        }
    }
}

async function fetchAndPatchJSON(url) {
    let resp = await fetch(url);
    const config = await resp.json();
    rewriteBaseUrlRecursive(config);
    return config;
}

async function loadData() {
    const datasources = await fetchAndPatchJSON(`${dir}/datasources.json`);
    const config = await fetchAndPatchJSON(`${dir}/state.json`);
    const views = await fetchAndPatchJSON(`${dir}/views.json`);
    const dataLoader = {
        function: getLocalCompressedBinaryDataLoader(datasources, dir),
        viewLoader: async (view) => views[view]
    }
    const cm = new ChartManager("app1", datasources, dataLoader, config);

};