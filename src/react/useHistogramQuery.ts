import type { HistogramInput, HistogramMessage } from '@/datastore/rawHistogramWorker';
import HistogramWorker from '@/datastore/rawHistogramWorker?worker';
import { useQuery } from '@tanstack/react-query';
import type { TypedArray } from 'd3';
import { useId } from 'react';

export function createHistogramWorker({data, min, max, bins, isInt32}: HistogramInput) {
    return new Promise<number[]>((resolve, reject) => {
        try {
            const worker = new HistogramWorker();

            worker.onmessage = (event) => {
                resolve(event.data);
                worker.terminate();
            };

            worker.onerror = (error) => {
                reject(error);
                worker.terminate();
            };

            const sharedBuffer = new SharedArrayBuffer(data.length * 4);
            new Float32Array(sharedBuffer).set(data);

            worker.postMessage({
                data: sharedBuffer,
                min,
                max,
                bins,
                isInt32,
            } satisfies HistogramMessage);
        } catch (err) {
            reject(err);
        }
    });
}

export function useHistogramQuery(
    //! todo: different data types, log scale...
    data: TypedArray,
    domain: [number, number],
    bins = 100,
    // should we set this with IntersectionObserver related state?
    enabled = true,
    name = "",
) {
    const id = useId();
    return useQuery({
        queryKey: [`histogram '${name}'`, id, domain[0], domain[1], bins, data.length],
        queryFn: () => createHistogramWorker({data, min: domain[0], max: domain[1], bins, isInt32: false}),
        enabled: enabled && data.length > 0,
        staleTime: Number.POSITIVE_INFINITY, // Data won't get stale unless explicitly invalidated
        // Only refetch when data or domain actually changes content
        structuralSharing: false,
    });
}
