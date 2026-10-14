import { type HistogramScaleMode, type NumericRange, resolveAutoHistogramXScaleFromValues } from "@/lib/utils";

/**
 * Upper bound on the number of raster values histogrammed per channel. Stats rasters come
 * from the coarsest pyramid level, which for an unpyramided image is the full-resolution
 * plane — tens of millions of pixels. A million samples is far more than 200 bins need.
 */
export const MAX_RASTER_HISTOGRAM_SAMPLES = 1 << 20;

export type RasterHistogramRequest = {
    values: Float64Array;
    domain: NumericRange;
    bins: number;
    xScaleMode: "auto" | HistogramScaleMode;
};

export type RasterHistogram = {
    counts: number[];
    edges: number[];
    xScale: HistogramScaleMode;
};

/**
 * Copy `values` into a transferable buffer, subsampling to at most `maxSamples`.
 *
 * Each sample is taken from a jittered position within its stride window rather than at a
 * fixed offset, so a stride that divides the image width can't alias onto a few columns.
 * The jitter is seeded, so the same raster always yields the same histogram.
 */
export function sampleRasterValues(values: ArrayLike<number>, maxSamples = MAX_RASTER_HISTOGRAM_SAMPLES) {
    const n = values.length;
    if (n <= maxSamples) {
        const out = new Float64Array(n);
        for (let i = 0; i < n; i++) out[i] = values[i];
        return out;
    }
    const out = new Float64Array(maxSamples);
    const stride = n / maxSamples;
    let seed = 0x9e3779b9;
    for (let i = 0; i < maxSamples; i++) {
        // xorshift32
        seed ^= seed << 13;
        seed ^= seed >>> 17;
        seed ^= seed << 5;
        const jitter = ((seed >>> 0) / 0x100000000) * stride;
        out[i] = values[Math.min(n - 1, Math.floor(i * stride + jitter))];
    }
    return out;
}

// d3.scaleSymlog with its default constant of 1.
const symlog = (x: number) => Math.sign(x) * Math.log1p(Math.abs(x));
const symexp = (y: number) => Math.sign(y) * Math.expm1(Math.abs(y));

/**
 * Count `values` into `bins` equal-width bins over `domain` — equal in symlog space when
 * `xScale` is "log". Non-finite and out-of-domain values are ignored. Always returns
 * exactly `bins` counts and `bins + 1` edges.
 */
export function binRasterValues(
    values: ArrayLike<number>,
    domain: NumericRange,
    bins: number,
    xScale: HistogramScaleMode,
): { counts: number[]; edges: number[] } {
    const min = domain[0];
    const max = domain[0] === domain[1] ? domain[0] + 1 : domain[1];
    const log = xScale === "log";
    const lo = log ? symlog(min) : min;
    const span = (log ? symlog(max) : max) - lo;
    const counts = new Array<number>(bins).fill(0);
    for (let i = 0; i < values.length; i++) {
        const value = values[i];
        if (!(value >= min && value <= max)) continue;
        const t = ((log ? symlog(value) : value) - lo) / span;
        counts[Math.min(bins - 1, Math.floor(t * bins))]++;
    }
    const edges = Array.from({ length: bins + 1 }, (_, i) => {
        const position = lo + (span * i) / bins;
        return log ? symexp(position) : position;
    });
    return { counts, edges };
}

export function computeRasterHistogram({ values, domain, bins, xScaleMode }: RasterHistogramRequest): RasterHistogram {
    const xScale = xScaleMode === "auto" ? resolveAutoHistogramXScaleFromValues(domain, values) : xScaleMode;
    return { ...binRasterValues(values, domain, bins, xScale), xScale };
}

/**
 * Compute a raster histogram in a worker. The raster is sampled on the calling thread (bounded
 * by `MAX_RASTER_HISTOGRAM_SAMPLES` reads) and the sample buffer is transferred, not copied.
 */
export function queryRasterHistogram(
    values: ArrayLike<number>,
    request: Omit<RasterHistogramRequest, "values">,
    signal?: AbortSignal,
) {
    if (values.length === 0) {
        return Promise.resolve(computeRasterHistogram({ ...request, values: new Float64Array() }));
    }
    return new Promise<RasterHistogram>((resolve, reject) => {
        const worker = new Worker(new URL("./rasterHistogramWorker.ts", import.meta.url), { type: "module" });
        const finish = (callback: () => void) => {
            worker.terminate();
            signal?.removeEventListener("abort", onAbort);
            callback();
        };
        const onAbort = () => finish(() => reject(new DOMException("Histogram query aborted", "AbortError")));
        if (signal?.aborted) {
            onAbort();
            return;
        }
        signal?.addEventListener("abort", onAbort, { once: true });
        worker.onmessage = (event: MessageEvent<RasterHistogram>) => finish(() => resolve(event.data));
        worker.onerror = (event) =>
            finish(() => reject(event.error ?? new Error(event.message || "Raster histogram worker failed")));
        const sample = sampleRasterValues(values);
        const message: RasterHistogramRequest = { ...request, values: sample };
        worker.postMessage(message, [sample.buffer]);
    });
}
