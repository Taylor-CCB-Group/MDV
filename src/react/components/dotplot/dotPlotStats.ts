/** Category codes for each data row: `text` columns are Uint8, `text16` columns are Uint16. */
export type CategoryCodes = Uint8Array | Uint16Array;

export type DotPlotStatsRequest = {
    /** Row category code for each data row (e.g. cell type). */
    rowCodes: CategoryCodes;
    /** Number of row categories, i.e. the row column's `values.length`. */
    rowCount: number;
    /** One array of values per field (x axis), each with one value per data row. */
    fields: ArrayLike<number>[];
    /** Indices of the data rows to include, e.g. from `useFilteredIndices()`. Defaults to all rows. */
    rows?: ArrayLike<number>;
    /**
     * Splits rows into panels. `panels` lists the split category codes to show, in panel order,
     * and must not repeat a code. Data rows in other categories are ignored.
     * Without a split, all rows go to a single panel.
     */
    split?: {
        codes: CategoryCodes;
        panels: number[];
    };
    /** A value counts as expressed when it is greater than this. Defaults to 0. */
    threshold?: number;
};

export type DotPlotCell = {
    /** Sum of expressed values divided by `count`, so unexpressed rows count as zero. */
    mean: number;
    /** Percentage (0-100) of rows in the group with an expressed value. */
    frac: number;
};

export type DotPlotGroup = {
    /** Number of data rows in this row category (and panel) after filtering. */
    count: number;
    /** One cell per field, in `fields` order. */
    cells: DotPlotCell[];
};

export type DotPlotStats = {
    /** One entry per panel, each with one group per row category code. */
    panels: DotPlotGroup[][];
    /** Min and max `mean` over every non-empty group in every panel; null if all are empty. */
    meanRange: [number, number] | null;
    /** Largest `frac` over every group in every panel. */
    fracMax: number;
};

/**
 * Mean expression and percentage expressing for each (panel, row category, field), in one pass
 * over the data. Uses the same definitions as the legacy dot plot's `addSimpleMean`
 * (`catColWorker.js`), grouped additionally by the split column.
 *
 * Unlike the legacy version, `meanRange` ignores empty groups, which are not drawn.
 */
export function computeDotPlotStats({
    rowCodes,
    rowCount,
    fields,
    rows,
    split,
    threshold = 0,
}: DotPlotStatsRequest): DotPlotStats {
    const panelCount = split ? split.panels.length : 1;
    const fieldCount = fields.length;
    const groupCount = panelCount * rowCount;

    const splitCodes = split?.codes;
    // split category code -> panel index, -1 if the category isn't shown
    const panelOfCode = new Int32Array(split ? Math.max(-1, ...split.panels) + 1 : 0).fill(-1);
    split?.panels.forEach((code, panel) => {
        panelOfCode[code] = panel;
    });

    const counts = new Uint32Array(groupCount);
    const sums = new Float64Array(groupCount * fieldCount);
    const expressing = new Uint32Array(groupCount * fieldCount);

    const rowTotal = rows ? rows.length : rowCodes.length;
    for (let j = 0; j < rowTotal; j++) {
        const i = rows ? rows[j] : j;
        let panel = 0;
        if (splitCodes) {
            const code = splitCodes[i];
            panel = code < panelOfCode.length ? panelOfCode[code] : -1;
            if (panel < 0) continue;
        }
        const group = panel * rowCount + rowCodes[i];
        counts[group]++;
        const base = group * fieldCount;
        for (let f = 0; f < fieldCount; f++) {
            const v = fields[f][i];
            // false for NaN, so missing values count as unexpressed
            if (v > threshold) {
                sums[base + f] += v;
                expressing[base + f]++;
            }
        }
    }

    let meanMin = Number.POSITIVE_INFINITY;
    let meanMax = Number.NEGATIVE_INFINITY;
    let fracMax = 0;
    const panels: DotPlotGroup[][] = [];
    for (let p = 0; p < panelCount; p++) {
        const groups: DotPlotGroup[] = [];
        for (let r = 0; r < rowCount; r++) {
            const group = p * rowCount + r;
            const count = counts[group];
            const cells: DotPlotCell[] = [];
            for (let f = 0; f < fieldCount; f++) {
                const k = group * fieldCount + f;
                const mean = count === 0 ? 0 : sums[k] / count;
                const frac = count === 0 ? 0 : (expressing[k] / count) * 100;
                cells.push({ mean, frac });
                if (count > 0) {
                    meanMin = Math.min(meanMin, mean);
                    meanMax = Math.max(meanMax, mean);
                }
                fracMax = Math.max(fracMax, frac);
            }
            groups.push({ count, cells });
        }
        panels.push(groups);
    }

    return {
        panels,
        meanRange: meanMin <= meanMax ? [meanMin, meanMax] : null,
        fracMax,
    };
}

/**
 * Compute dot plot stats in a worker. Typed arrays backed by a SharedArrayBuffer (as MDV column
 * data is) are shared with the worker rather than copied.
 */
export function queryDotPlotStats(request: DotPlotStatsRequest, signal?: AbortSignal) {
    return new Promise<DotPlotStats>((resolve, reject) => {
        const worker = new Worker(new URL("./dotPlotStatsWorker.ts", import.meta.url), { type: "module" });
        const finish = (callback: () => void) => {
            worker.terminate();
            signal?.removeEventListener("abort", onAbort);
            callback();
        };
        const onAbort = () => finish(() => reject(new DOMException("Dot plot stats query aborted", "AbortError")));
        if (signal?.aborted) {
            onAbort();
            return;
        }
        signal?.addEventListener("abort", onAbort, { once: true });
        worker.onmessage = (event: MessageEvent<DotPlotStats>) => finish(() => resolve(event.data));
        worker.onerror = (event) =>
            finish(() => reject(event.error ?? new Error(event.message || "Dot plot stats worker failed")));
        worker.postMessage(request);
    });
}
