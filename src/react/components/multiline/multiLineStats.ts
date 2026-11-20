/** Category codes for each data row: `text` columns are Uint8, `text16` columns are Uint16. */
export type CategoryCodes = Uint8Array | Uint16Array;

export type MultiLineDensityRequest = {
    /** Category code for each data row: one curve per category. */
    categoryCodes: CategoryCodes;
    /** Number of categories, i.e. the category column's `values.length`. */
    categoryCount: number;
    /** The value of each data row. NaN values are left out. */
    values: ArrayLike<number>;
    /** Indices of the data rows to include, e.g. from `useFilteredIndices()`. Defaults to all rows. */
    rows?: ArrayLike<number>;
    /** Where to evaluate the density, in ascending order and in the units of `values`. */
    ticks: ArrayLike<number>;
    /** Half-width of the kernel, in the units of `values`. */
    bandwidth: number;
    /** Clamp values to this range first, as the legacy chart does when trimmed to a percentile. */
    clamp?: [number, number];
};

/** The density of one category at each tick. */
export type DensityCurve = {
    /** Category code. */
    category: number;
    density: Float64Array;
    /** Largest value in `density`: each curve is drawn scaled to its own peak. */
    max: number;
};

/** Curves in category order, leaving out categories with no density (no rows, or all NaN). */
export type MultiLineDensities = { curves: DensityCurve[] };

/** First index in ascending `ticks` whose value is at least `x`. */
function lowerBound(ticks: ArrayLike<number>, x: number) {
    let lo = 0;
    let hi = ticks.length;
    while (lo < hi) {
        const mid = (lo + hi) >>> 1;
        if (ticks[mid] < x) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

/**
 * Kernel density of each category's values (Epanechnikov kernel), as the legacy multi-line chart
 * computes it in `boxPlotWorker.js`: at each tick, the mean over the category's values of
 * `0.75 * (1 - u²) / bandwidth` where `u = (tick - value) / bandwidth` and `|u| <= 1`.
 *
 * Each value only adds to the ticks within `bandwidth` of it, so this is linear in the rows rather
 * than in rows × ticks. Contributions are added in row order, as the legacy worker does, so the
 * results match it exactly.
 */
export function computeMultiLineDensities({
    categoryCodes,
    categoryCount,
    values,
    rows,
    ticks,
    bandwidth,
    clamp,
}: MultiLineDensityRequest): MultiLineDensities {
    const tickCount = ticks.length;
    const sums = Array.from({ length: categoryCount }, () => new Float64Array(tickCount));
    const counts = new Float64Array(categoryCount);
    const k = bandwidth;
    const rowCount = rows ? rows.length : values.length;
    const [lo, hi] = clamp ?? [Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY];

    for (let r = 0; r < rowCount; r++) {
        const i = rows ? rows[r] : r;
        let v = values[i];
        if (Number.isNaN(v)) continue;
        const category = categoryCodes[i];
        if (category >= categoryCount) continue;
        v = v < lo ? lo : v > hi ? hi : v;
        counts[category]++;
        const sum = sums[category];
        // a tick just outside the kernel can still pass the test below by rounding, so start one early
        for (let j = Math.max(0, lowerBound(ticks, v - k) - 1); j < tickCount; j++) {
            const x = ticks[j];
            if (x - v > k) break;
            const u = (x - v) / k;
            if (Math.abs(u) <= 1) sum[j] += (0.75 * (1 - u * u)) / k;
        }
    }

    const curves: DensityCurve[] = [];
    for (let c = 0; c < categoryCount; c++) {
        const count = counts[c];
        if (count === 0) continue;
        const density = sums[c];
        let max = 0;
        for (let j = 0; j < tickCount; j++) {
            density[j] /= count;
            if (density[j] > max) max = density[j];
        }
        if (max === 0 || Number.isNaN(max)) continue;
        curves.push({ category: c, density, max });
    }
    return { curves };
}

/**
 * Compute densities in a worker. Each query has its own worker, so concurrent queries (one per
 * cell of a grid) each get their own answer. Typed arrays backed by a SharedArrayBuffer (as MDV
 * column data is) are shared with the worker rather than copied.
 */
export function queryMultiLineDensities(request: MultiLineDensityRequest, signal?: AbortSignal) {
    return new Promise<MultiLineDensities>((resolve, reject) => {
        const worker = new Worker(new URL("./multiLineStatsWorker.ts", import.meta.url), { type: "module" });
        const finish = (callback: () => void) => {
            worker.terminate();
            signal?.removeEventListener("abort", onAbort);
            callback();
        };
        const onAbort = () => finish(() => reject(new DOMException("Density query aborted", "AbortError")));
        if (signal?.aborted) {
            onAbort();
            return;
        }
        signal?.addEventListener("abort", onAbort, { once: true });
        worker.onmessage = (event: MessageEvent<MultiLineDensities>) => finish(() => resolve(event.data));
        worker.onerror = (event) =>
            finish(() => reject(event.error ?? new Error(event.message || "Density worker failed")));
        worker.postMessage(request);
    });
}
