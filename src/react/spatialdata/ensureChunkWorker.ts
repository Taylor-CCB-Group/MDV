import {
    enableWorkerChunkDecode,
    type EnableWorkerChunkDecodeOptions,
} from "zarrextra/workers";

let chunkWorkerEnabled = false;

/** Lazily enable zarrextra worker-pool chunk decode before spatial zarr reads. */
export function ensureChunkWorker(options?: EnableWorkerChunkDecodeOptions) {
    if (chunkWorkerEnabled) return;
    enableWorkerChunkDecode(options);
    chunkWorkerEnabled = true;
}
