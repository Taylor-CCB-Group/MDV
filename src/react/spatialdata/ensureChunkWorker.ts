import {
    enableWorkerChunkDecode,
    type EnableWorkerChunkDecodeOptions,
} from "zarrextra/workers";

let chunkWorkerEnabled = false;

/** Lazily enable zarrextra worker-pool chunk decode before spatial zarr reads.
 * This is currently redundant but in future we may revive to pass in options
 */
export function ensureChunkWorker(options?: EnableWorkerChunkDecodeOptions) {
    if (chunkWorkerEnabled) return;
    // we should get rid of this and let spatialdata enable by default.
    // enableWorkerChunkDecode(options);
    console.log("redundant call to enableWorkerChunkDecode which should now be done by default")
    chunkWorkerEnabled = true;
}
