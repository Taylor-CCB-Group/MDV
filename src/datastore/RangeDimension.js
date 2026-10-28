import Dimension from "./Dimension.js";
import { RangeFilterClient } from "./rangeFilterClient";
import { evaluateFilterSpec } from "./rangeFilterEval";

class RangeDimension extends Dimension {
    /**
     * @param {DataStore} parent
     */
    constructor(parent) {
        super(parent);
        this.worker = new Worker(
            new URL("./binWorker.js?v=1", import.meta.url),
        );
    }

    filterSquare(args, columns) {
        return this._commitSpec({
            kind: "square",
            data1: this._columnData(columns[0]),
            data2: this.parent.columnIndex[columns[1]].data,
            // plain copies: args may be observable proxies, which can't be posted to a worker
            range1: [args.range1[0], args.range1[1]],
            range2: [args.range2[0], args.range2[1]],
        });
    }

    /**
     * Polygon containment filter on x/y columns.
     * 
     * args can be either:
     *   - an array of polygon vertices
     *   - { points, regionField, regionValue } to additionally restrict rows to a
     *     specific region value.
     * 
     * @param {Array<[number, number]> | { points: Array<[number, number]>, regionField: string, regionValue: string }} args 
     *
     */
    filterPoly(args, columns) {
        const points = Array.isArray(args) ? args : args.points;
        const regionField = Array.isArray(args) ? undefined : args.regionField;
        const regionValue = Array.isArray(args) ? undefined : args.regionValue;
        /** @type {import("./rangeFilterEval").PolyFilterSpec} */
        const spec = {
            kind: "poly",
            data1: this._columnData(columns[0]),
            data2: this.parent.columnIndex[columns[1]].data,
            points: points.map((pt) => [pt[0], pt[1]]),
        };
        // for region-scoped gates, resolve the wanted region value's index once, up front
        if (regionField && regionValue != null) {
            const regionCol = this.parent.columnIndex[regionField];
            spec.region = {
                data: regionCol?.data ?? null,
                index: regionCol?.values ? regionCol.values.indexOf(regionValue) : -1,
            };
        }
        return this._commitSpec(spec);
    }

    filterRange(args, columns) {
        return this._commitSpec({
            kind: "range",
            data: this.parent.columnIndex[columns[0]].data,
            min: args.min,
            max: args.max,
        });
    }

    /** filterSquare and filterPoly accept a typed array in place of the first column name */
    _columnData(column) {
        return typeof column === "string" ? this.parent.columnIndex[column].data : column;
    }

    /**
     * Evaluates a filter spec on this dimension's filter worker, or synchronously below
     * {@link RangeDimension.asyncFilterMinRows}.
     * @param {import("./rangeFilterEval").FilterSpec} spec
     */
    _commitSpec(spec) {
        const size = this.parent.size;
        if (size < RangeDimension.asyncFilterMinRows) {
            const mask = new Uint8Array(size);
            evaluateFilterSpec(spec, size, mask);
            return this.commitFilter(mask);
        }
        this._filterClient ??= new RangeFilterClient(
            () => new Worker(new URL("./rangeFilterWorker.ts", import.meta.url), { type: "module" }),
        );
        return this.commitFilter(this._filterClient.evaluate(spec, size));
    }

    cancelPendingFilter() {
        super.cancelPendingFilter();
        this._filterClient?.cancelQueued();
    }

    getBins(callback, column, config = {}) {
        const col = this.parent.columnIndex[column];
        config.bins = config.bins === undefined ? 10 : config.bins;
        const [min, max] = this.parent.getMinMaxForColumn(column);
        config.min = config.min === undefined ? min : config.min;
        config.max = config.max === undefined ? max : config.max;

        const t = performance.now();
        const action = (e) => {
            console.log(
                `calculate bins ${col.name} : ${performance.now() - t}`,
            );
            callback(e.data);
        };
        this.worker.onmessage = action;

        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            [col.buffer, col.datatype],
            config,
        ]);
    }
    getBinsAsync(column, config = {}) {
        // wasted hours on this, gave up and made a simpler version for unfiltered histogram.
        return new Promise((resolve) => {
            this.getBins(resolve, column, config);
        });
    }

    destroy() {
        super.destroy();
        this.worker.terminate();
        this._filterClient?.terminate();
    }
}

/**
 * Datasets with fewer rows than this filter synchronously, so the result is visible before
 * `filter()` returns. Larger ones evaluate on the dimension's filter worker; the main thread then
 * only commits the mask.
 *
 * The worker round trip itself is cheap (~0.1ms main thread, ~0.3ms extra latency), so this is
 * about keeping small datasets synchronous rather than about overhead: at 100k rows a sync range
 * filter takes <1ms and a 24-vertex polygon ~7ms; at 1M those are ~5ms and ~60ms. Set to 0 to
 * always use the worker (tests), or Infinity to never use it.
 */
RangeDimension.asyncFilterMinRows = 100_000;

Dimension.types["range_dimension"] = RangeDimension;

export default RangeDimension;
