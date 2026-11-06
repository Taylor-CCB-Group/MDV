import { Matrix4, type NumericArray } from "@math.gl/core";

/**
 * The frame spatial statistics are computed in, and how MDV's position columns get there.
 *
 * Statistics take lengths (rMax, bin width), so the cells and the observation window must be in the
 * same physical frame — µm. `toAnalysis` maps the data-column frame (what `regions.position_fields`,
 * gate geometry and the live lasso selection are all stored in) into it. It is a 4×4 so a z column
 * and a differing z spacing slot in later; the 2-D statistics use its x/y rows only.
 *
 * Domains and cells must both go through the same frame: build one per computation and pass it to
 * both `DomainSource.toDomain` and the cell gather.
 */
export type AnalysisFrame = {
    /** Data columns → analysis frame, column-major (math.gl `Matrix4` layout). */
    readonly toAnalysis: readonly number[];
    /** Unit of the analysis frame, verbatim, for display. */
    readonly unit: string;
    /** Where the transform came from, for debug output and to make stopgaps visible. */
    readonly provenance: string;
};

export function frameFromMatrix(m: Matrix4 | NumericArray, unit: string, provenance: string): AnalysisFrame {
    return { toAnalysis: Array.from(m), unit, provenance };
}

export function identityFrame(unit = "data units"): AnalysisFrame {
    return frameFromMatrix(new Matrix4().identity(), unit, "identity");
}

/**
 * `frame` followed by a stated physical size of its unit — for stores that leave the unit unstated
 * (SpatialData's "unit" placeholder; Xenium writes it everywhere). The caller says where the number
 * came from; it is carried into `provenance` so a result never silently claims µm.
 */
export function withPhysicalUnit(frame: AnalysisFrame, micrometresPerUnit: number, why: string): AnalysisFrame {
    const s = micrometresPerUnit;
    const m = new Matrix4().scale([s, s, s]).multiplyRight(new Matrix4(frame.toAnalysis));
    return frameFromMatrix(m, "µm", `${frame.provenance}, ×${s} µm/unit (${why})`);
}

/** The x/y part of `toAnalysis` as `[a, b, c, d, tx, ty]`: `x' = a·x + c·y + tx`, `y' = b·x + d·y + ty`. */
export type Affine2 = readonly [number, number, number, number, number, number];

export function affine2(frame: AnalysisFrame): Affine2 {
    const m = frame.toAnalysis;
    return [m[0], m[1], m[4], m[5], m[12], m[13]];
}

export function applyAffine2([a, b, c, d, tx, ty]: Affine2, x: number, y: number): [number, number] {
    return [a * x + c * y + tx, b * x + d * y + ty];
}

/**
 * Today's only reachable source of physical scale for a datastore: `regions.scale` /
 * `regions.scale_unit`, one isotropic factor for every region (see `MDVProject.set_region_scale`).
 *
 * This is the stopgap the cross-PCF slice started from, now behind the frame seam. Its successor is
 * the SpatialData coordinate-system transform of the element the table annotates — see
 * `docs/` / the session notes for what is and isn't reachable from the datastore yet.
 */
export function regionsScaleFrame(regions: { scale?: number; scale_unit?: string } | undefined): AnalysisFrame {
    const s = regions?.scale;
    if (!s || !Number.isFinite(s)) return identityFrame();
    return frameFromMatrix(new Matrix4().scale([s, s, s]), regions?.scale_unit ?? "µm", "regions.scale");
}
