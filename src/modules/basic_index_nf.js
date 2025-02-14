import "microtip/microtip.min.css";
import "nouislider/dist/nouislider.min.css";
import "../utilities/css/ContextMenu.css";
import "../charts/css/charts.css";
import "../webgl/css/wgl2di.css";
import "../table/css/slickgrid.css";
import "../browser/css/browser.css";
import ChartManager from "../charts/ChartManager.js";
import {
    getArrayBufferDataLoader,
    getLocalCompressedBinaryDataLoader,
    processArrayBuffer,
    decompressData,
} from "../dataloaders/DataLoaders";
window.mdv = {
    ChartManager,
    getArrayBufferDataLoader,
    getLocalCompressedBinaryDataLoader,
    processArrayBuffer,
    decompressData,
};
