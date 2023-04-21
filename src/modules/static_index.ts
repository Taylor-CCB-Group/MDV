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
const dir = urlParams.get('dir') || '/data/ytrap';//http://localhost:8082/';

async function loadData() {
    let resp = await fetch(`${dir}/datasources.json`);
    const datasources = await resp.json();
    resp = await fetch(`${dir}/state.json`);
    const config = await resp.json();
    resp = await fetch(`${dir}/views.json`);
    const views = await resp.json();
    const dataLoader = {
        function: getLocalCompressedBinaryDataLoader(datasources, dir),
        viewLoader: async (view) => views[view]
    }
    const cm = new ChartManager("app1", datasources, dataLoader, config);

};