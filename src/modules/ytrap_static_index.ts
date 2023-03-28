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

const ytrap = 'http://localhost:8080/';

async function loadData() {
    let resp = await fetch(`${ytrap}/datasources.json`);
    const datasources = await resp.json();
    resp = await fetch(`${ytrap}/state.json`);
    const config = await resp.json();
    resp = await fetch(`${ytrap}/views.json`);
    const views = await resp.json();
    const dataLoader = {
        function: getLocalCompressedBinaryDataLoader(datasources, ytrap),
        viewLoader: async (view) => views[view]
    }
    const cm = new ChartManager("app1", datasources, dataLoader, config);

};