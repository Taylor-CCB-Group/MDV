import type { NumericArrayType, NumericColumnData } from "@/lib/columnTypeHelpers";

export type HistogramWorkerInput = {
    data: SharedArrayBuffer;
    min: number;
    max: number;
    bins: number;
    arrayType: NumericArrayType;
    byteOffset: number;
    length: number;
};

export function createHistogram(
    originalData: NumericColumnData,
    min: number,
    max: number,
    bins: number,
    indices?: ArrayLike<number>,
) {
    if (bins <= 0 || max <= min) return new Array(bins).fill(0);
    const histogram = new Array(bins).fill(0);
    const binWidth = (max - min) / bins;
    const addValue = (value: number) => {
        if (Number.isNaN(value)) return;
        const bin = value === max ? bins - 1 : Math.floor((value - min) / binWidth);
        if (bin >= 0 && bin < bins) {
            histogram[bin]++;
        }
    };

    if (indices) {
        for (let i = 0; i < indices.length; i++) {
            addValue(originalData[indices[i]]);
        }
        return histogram;
    }

    for (let i = 0; i < originalData.length; i++) {
        addValue(originalData[i]);
    }
    return histogram;
}

export function queryHistogramWorker(input: HistogramWorkerInput, signal?: AbortSignal) {
    return new Promise<number[]>((resolve, reject) => {
        const worker = new Worker(new URL("../../datastore/rawHistogramWorker.ts", import.meta.url));
        let settled = false;

        const finish = (callback: () => void) => {
            if (settled) return;
            settled = true;
            worker.terminate();
            callback();
        };

        const onAbort = () => {
            finish(() => reject(new DOMException("Histogram query aborted", "AbortError")));
        };

        worker.onmessage = (event) => {
            finish(() => resolve(event.data));
        };

        worker.onerror = (event) => {
            finish(() => reject(event.error ?? new Error(event.message || "Histogram worker failed")));
        };

        signal?.addEventListener("abort", onAbort, { once: true });
        if (signal?.aborted) {
            onAbort();
            return;
        }

        worker.postMessage(input);
    });
}
