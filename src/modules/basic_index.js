
import "./all_css"
import ChartManager from '../charts/ChartManager.js';
import { getArrayBufferDataLoader,getLocalCompressedBinaryDataLoader,processArrayBuffer} from "../dataloaders/DataLoaders.js";
window.mdv={
    ChartManager,
    getArrayBufferDataLoader,
    getLocalCompressedBinaryDataLoader,
    processArrayBuffer
}