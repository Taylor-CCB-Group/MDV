import "./all_css";
import ChartManager from "../charts/ChartManager.js";
import {
    decompressData,
    getArrayBufferDataLoader,
    getLocalCompressedBinaryDataLoader,
    processArrayBuffer,
} from "../dataloaders/DataLoaders";

window.mdv = {
    ChartManager,
    getArrayBufferDataLoader,
    getLocalCompressedBinaryDataLoader,
    processArrayBuffer,
    decompressData,
};
