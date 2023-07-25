import "microtip/microtip.min.css";
import 'nouislider/dist/nouislider.min.css'
import "../utilities/css/ContextMenu.css";
import "../charts/css/charts.css";
import "../css/fontawesome-5.15.3/all.min.css";  //remove
import "../webgl/css/wgl2di.css";
import "../table/css/slickgrid.css";
import "../charts/CustomBoxPlot.js";
import "../charts/MultiBoxPlot.js";
import "../charts/ImageTableChart.js";
import "../charts/SingleHeatMap.js";
import "../charts/MultiPieChart.js";
import "../charts/SingleSeriesChart.js";
import "../charts/ImageBox.js";
import "../charts/SelectionDialog";
import "../charts/SingleHeatMapGroup.js";
import "../charts/VivScatterPlot.js";
import "../charts/CellNetworkChart.js";
import "../charts/CellRadialChart.js";
import "../charts/RowSummaryBox.js";
import "../charts/ImageScatterChart.js";
import "../charts/WordCloudChart.js";


import "../browser/css/browser.css";
import "../charts/GenomeBrowser.js";
import { createEl } from "../utilities/Elements";
import connectWebsocket from '../utilities/WebSocket';

//objects to attach to the window object
import ChartManager from '../charts/ChartManager.js';
import FileUploadDialog from "../utilities/FileUploadDialog.js";
import { processArrayBuffer, getArrayBufferDataLoader,getLocalCompressedBinaryDataLoader} from "../dataloaders/DataLoaders.js";
window.ChartManager = ChartManager;
window.scFileUploadDialog= FileUploadDialog;
window.processArrayBuffer = processArrayBuffer;
window.getArrayBufferDataLoader= getArrayBufferDataLoader;
window.createEl= createEl;
window.getLocalCompressedBinaryDataLoader = getLocalCompressedBinaryDataLoader;
//window.cm will be set in link.js, at least for now
setTimeout(()=>connectWebsocket('http://localhost:5050', window.cm), 1000);
