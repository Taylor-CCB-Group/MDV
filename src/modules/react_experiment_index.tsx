import ReactTest from "./react_test";
console.log(ReactTest);
import "microtip/microtip.min.css";
import "../css/fontawesome-5.15.3/all.min.css";
import 'nouislider/dist/nouislider.min.css'
import "../utilities/css/ContextMenu.css";
import "../charts/css/charts.css";
import "../webgl/css/wgl2di.css";
import "../table/css/slickgrid.css";
// import ChartManager from '../charts/ChartManager.js';
import "../charts/VivScatterPlot";
import "../webgl/VivVolume.js";
import chartManager from "../../examples/obvios_example";
import AnnotationDialog from "../charts/dialogs/AnnotationDialog";
document.title = "MDV - OBVioS";
(window as any).chartManager = chartManager;

const ds = chartManager.dsIndex['test'];
chartManager.addMenuIcon(ds.name, "fas fa-tags", "Tag annotation", () => { new AnnotationDialog(ds.dataStore) });