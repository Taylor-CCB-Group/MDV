import Dimension from "./Dimension.js";

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
        let data1 = null;
        if (typeof columns[0] !== "string") {
            data1 = columns[0];
        } else {
            data1 = this.parent.columnIndex[columns[0]].data;
        }

        const data2 = this.parent.columnIndex[columns[1]].data;
        const range1 = args.range1;
        const range2 = args.range2;
        const predicate = (i) => {
            const v1 = data1[i];
            const v2 = data2[i];
            return (
                v1 >= range1[0] &&
                v1 <= range1[1] &&
                v2 >= range2[0] &&
                v2 <= range2[1] &&
                Number.isFinite(v1) &&
                Number.isFinite(v2)
            );
        };
        return this.filterPredicate({ predicate }, columns);
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

        // Compute bounding box
        let minX = Number.MAX_VALUE;
        let minY = Number.MAX_VALUE;
        let maxX = Number.MIN_VALUE;
        let maxY = Number.MIN_VALUE;
        for (const pt of points) {
            minX = Math.min(minX, pt[0]);
            maxX = Math.max(maxX, pt[0]);
            minY = Math.min(minY, pt[1]);
            maxY = Math.max(maxY, pt[1]);
        }
        let data1 = null;
        if (typeof columns[0] !== "string") {
            data1 = columns[0];
        } else {
            data1 = this.parent.columnIndex[columns[0]].data;
        }
        const data2 = this.parent.columnIndex[columns[1]].data;

        // For region-scoped gates, resolve the region column and the index of the
        // wanted region value once — before the per-row loop — so it's not repeated.
        const isRegionScoped = Boolean(regionField && regionValue != null);
        const regionCol = isRegionScoped ? this.parent.columnIndex[regionField] : null;
        const regionData = regionCol?.data ?? null;
        const wantedRegionIndex = regionCol?.values ? regionCol.values.indexOf(regionValue) : -1;

        const predicate = (i) => {
            const x = data1[i];
            const y = data2[i];
            if (
                x < minX ||
                x > maxX ||
                y < minY ||
                y > maxY ||
                !Number.isFinite(x) ||
                !Number.isFinite(y)
            ) {
                return false;
            }

            // For region-scoped gates, check that this row belongs to the gate's region
            if (isRegionScoped) {
                if (!regionData || wantedRegionIndex === -1) return false;
                if (regionData[i] !== wantedRegionIndex) return false;
            }

            // Ray-casting algorithm: count edge crossings to determine if point is inside polygon.
            let inside = false;
            for (let curr = 0, prev = points.length - 1; curr < points.length; prev = curr++) {
                const xi = points[curr][0];
                const yi = points[curr][1];
                const xj = points[prev][0];
                const yj = points[prev][1];
                const intersect =
                    yi > y !== yj > y &&
                    x < ((xj - xi) * (y - yi)) / (yj - yi) + xi;
                if (intersect) inside = !inside;
            }

            return inside;
        };

        return this.filterPredicate({ predicate }, columns);
    }

    filterRange(args, columns) {
        const min = args.min;
        const max = args.max;
        const arr = this.parent.columnIndex[columns[0]].data;
        // performance seems similar to non-predicate version
        const predicate = (i) => {
            const v = arr[i];
            return v >= min && v <= max && !Number.isNaN(v);
        };
        return this.filterPredicate({ predicate }, columns);
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
    }
}

Dimension.types["range_dimension"] = RangeDimension;

export default RangeDimension;
