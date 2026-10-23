/**
 * Pure, serialisable evaluation of RangeDimension's range / square / polygon filters.
 *
 * Shared by `rangeFilterWorker.ts` (off-main-thread evaluation) and RangeDimension's sync
 * path, so both produce identical results: inclusive bounds, NaN (and, for square and
 * polygon, any non-finite value) excluded, ray-casting for polygons, optional region scoping.
 *
 * Column arrays are passed as typed arrays; when they are views over a SharedArrayBuffer (as
 * DataStore columns are) posting a spec to a worker shares the memory rather than copying it.
 */

/** Numeric column values as stored by DataStore (or a caller-supplied typed array). */
export type NumericArray = ArrayLike<number>;

export type RangeFilterSpec = {
    kind: "range";
    data: NumericArray;
    min: number;
    max: number;
};

export type SquareFilterSpec = {
    kind: "square";
    data1: NumericArray;
    data2: NumericArray;
    range1: readonly [number, number];
    range2: readonly [number, number];
};

export type PolyFilterSpec = {
    kind: "poly";
    data1: NumericArray;
    data2: NumericArray;
    points: ReadonlyArray<readonly [number, number]>;
    /**
     * Restricts rows to one region. `data` null or `index` -1 (column or value not found)
     * excludes every row, matching the sync filter.
     */
    region?: { data: NumericArray | null; index: number };
};

export type FilterSpec = RangeFilterSpec | SquareFilterSpec | PolyFilterSpec;

export type FilterRequest = {
    id: number;
    spec: FilterSpec;
    size: number;
    /** receives the include mask: 1 keeps the row */
    output: SharedArrayBuffer;
};

export type FilterResponse = { id: number; error?: string };

/**
 * Writes the include mask (1 keeps the row, 0 drops it) for rows `[0, size)` into `out`.
 *
 * Each kind has its own loop rather than sharing a per-row callback: a callback site that sees
 * several kinds (as the worker's does) stops being inlined, which made worker polygons ~2x slower.
 */
export function evaluateFilterSpec(spec: FilterSpec, size: number, out: Uint8Array) {
    if (out.length < size) {
        throw new Error(`filter output has ${out.length} rows, expected ${size}`);
    }
    switch (spec.kind) {
        case "range":
            return evaluateRange(spec, size, out);
        case "square":
            return evaluateSquare(spec, size, out);
        case "poly":
            return evaluatePoly(spec, size, out);
    }
}

/** Worker message handler, also driven in-process by tests. */
export function handleFilterRequest(request: FilterRequest): FilterResponse {
    try {
        evaluateFilterSpec(request.spec, request.size, new Uint8Array(request.output));
        return { id: request.id };
    } catch (e) {
        return { id: request.id, error: e instanceof Error ? e.message : String(e) };
    }
}

function evaluateRange({ data, min, max }: RangeFilterSpec, size: number, out: Uint8Array) {
    for (let i = 0; i < size; i++) {
        const v = data[i];
        // inclusive; NaN fails both comparisons
        out[i] = v >= min && v <= max ? 1 : 0;
    }
}

function evaluateSquare({ data1, data2, range1, range2 }: SquareFilterSpec, size: number, out: Uint8Array) {
    const [min1, max1] = range1;
    const [min2, max2] = range2;
    for (let i = 0; i < size; i++) {
        const v1 = data1[i];
        const v2 = data2[i];
        // inclusive; infinities are excluded even when a bound is infinite
        out[i] =
            v1 >= min1 && v1 <= max1 && v2 >= min2 && v2 <= max2 && Number.isFinite(v1) && Number.isFinite(v2) ? 1 : 0;
    }
}

function evaluatePoly({ data1, data2, points, region }: PolyFilterSpec, size: number, out: Uint8Array) {
    if (region && (!region.data || region.index === -1)) {
        out.fill(0, 0, size);
        return;
    }
    const regionData = region?.data;
    const regionIndex = region?.index;
    const n = points.length;
    const px = Float64Array.from(points, (pt) => pt[0]);
    const py = Float64Array.from(points, (pt) => pt[1]);
    // bounding box as a cheap pre-test; Number.MIN_VALUE (smallest positive) for the max seeds
    // is historical and only ever widens the box, so ray-casting still decides
    let minX = Number.MAX_VALUE;
    let minY = Number.MAX_VALUE;
    let maxX = Number.MIN_VALUE;
    let maxY = Number.MIN_VALUE;
    for (let k = 0; k < n; k++) {
        minX = Math.min(minX, px[k]);
        maxX = Math.max(maxX, px[k]);
        minY = Math.min(minY, py[k]);
        maxY = Math.max(maxY, py[k]);
    }
    for (let i = 0; i < size; i++) {
        const x = data1[i];
        const y = data2[i];
        if (x < minX || x > maxX || y < minY || y > maxY || !Number.isFinite(x) || !Number.isFinite(y)) {
            out[i] = 0;
            continue;
        }
        if (regionData && regionData[i] !== regionIndex) {
            out[i] = 0;
            continue;
        }
        // ray-casting: count edge crossings to decide whether the point is inside
        let inside = false;
        for (let curr = 0, prev = n - 1; curr < n; prev = curr++) {
            const xi = px[curr];
            const yi = py[curr];
            const xj = px[prev];
            const yj = py[prev];
            if (yi > y !== yj > y && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
        }
        out[i] = inside ? 1 : 0;
    }
}
