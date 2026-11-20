import type { Domain } from "intraspatial/spatial/domain";
import { crossPCFMatrixBinned, type LabelledCells, type PcfMatrixBinnedResult } from "intraspatial/spatial/pcf";
import { type AnalysisFrame, affine2 } from "./analysisFrame";

/**
 * Column arrays exactly as DataStore holds them (`column.data`): a view over its (Shared)ArrayBuffer,
 * so posting one shares the memory rather than copying it when the buffer is shared.
 *
 * The typed array is passed rather than a buffer + `datatype`, because the datatype → array mapping
 * is DataStore's to decide and isn't the obvious one: `integer` is held as `Float32Array` today
 * (see `DataStructureTypes`), and `datatypes.ts` still says `Int32Array`. Reading the array DataStore
 * built means this follows whatever that mapping becomes.
 */
export type NumericArray = Float32Array | Int32Array;
/** Category codes: `text` → Uint8, `text16` → Uint16. */
export type CodeArray = Uint8Array | Uint16Array;

/** What the cross-PCF worker receives. Plain data, so it posts as-is. */
export type CrossPcfRequest = {
    x: NumericArray;
    y: NumericArray;
    /** Categorical column: its codes are the type ids. */
    type: { codes: CodeArray; nTypes: number };
    /** Keep only rows whose region code equals `code`. */
    region?: { codes: CodeArray; code: number };
    /** Data columns → analysis frame; `domain` is already in the analysis frame. */
    frame: AnalysisFrame;
    domain: Domain;
    rMax: number;
    nBins: number;
    edgeCorrected: boolean;
};

export type CrossPcfResponse = PcfMatrixBinnedResult & {
    /** Rows in the domain's bounding box handed to intraspatial (before its exact inside test). */
    candidates: number;
    ms: number;
};

/**
 * The cells handed to intraspatial: positions already in the analysis frame, as typed arrays.
 *
 * intraspatial's `LabelledCells` is `readonly number[]` today, so this is the one place MDV copies to
 * a JS array. Once it takes `Coords` / integer typed arrays, this becomes `({ xs, ys, typeId })`.
 */
export function labelledCells(xs: Float64Array, ys: Float64Array, typeId: Uint16Array): LabelledCells {
    return { xs: Array.from(xs), ys: Array.from(ys), typeId: Array.from(typeId) };
}

/**
 * Gather the rows that can be inside `domain`, mapped through `frame`. The bbox test is a cheap
 * prefilter; `crossPCFMatrixBinned` applies the exact inside test.
 *
 * The transformed positions have to be materialised somewhere (the transform is not the identity),
 * so they go straight into typed arrays sized for the worst case and are trimmed with `subarray`.
 */
export function gatherCells(req: Pick<CrossPcfRequest, "x" | "y" | "type" | "region" | "frame" | "domain">) {
    const { x, y } = req;
    const type = req.type.codes;
    const region = req.region?.codes;
    const regionCode = req.region?.code;
    const [minX, minY, maxX, maxY] = req.domain.bbox;
    const [a, b, c, d, tx, ty] = affine2(req.frame);
    const nTypes = req.type.nTypes;

    const n = x.length;
    const xs = new Float64Array(n);
    const ys = new Float64Array(n);
    const typeId = new Uint16Array(n);
    let k = 0;
    for (let i = 0; i < n; i++) {
        if (region && region[i] !== regionCode) continue;
        const px = a * x[i] + c * y[i] + tx;
        const py = b * x[i] + d * y[i] + ty;
        if (!(px >= minX && px <= maxX && py >= minY && py <= maxY)) continue; // also drops NaN
        const t = type[i];
        if (t >= nTypes) continue; // missing-value code beyond the values list
        xs[k] = px;
        ys[k] = py;
        typeId[k] = t;
        k++;
    }
    return { xs: xs.subarray(0, k), ys: ys.subarray(0, k), typeId: typeId.subarray(0, k) };
}

/** Cross-PCF matrix of the cells inside `domain`, with `domain` as the observation window. */
export function crossPcfInDomain(req: CrossPcfRequest): CrossPcfResponse {
    const t0 = performance.now();
    const { xs, ys, typeId } = gatherCells(req);
    const result = crossPCFMatrixBinned(labelledCells(xs, ys, typeId), {
        bbox: req.domain.bbox, // required by the type even with `domain`; ignored then
        domain: req.domain,
        rMax: req.rMax,
        nBins: req.nBins,
        nTypes: req.type.nTypes,
        edgeCorrected: req.edgeCorrected,
    });
    return { ...result, candidates: xs.length, ms: performance.now() - t0 };
}
