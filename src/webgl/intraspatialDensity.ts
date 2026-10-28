import type { ContourLayerProps } from "@/react/contour_state";
import type { DensityLayerProps } from "@intraspatial/deck";
import type { Matrix4 } from "@math.gl/core";

/**
 * Map MDV's heatmap contour props onto @intraspatial/deck's `DensityLayer`, which draws
 * contours on WebGPU (where `TriangleLayerContours`, being GLSL, cannot).
 *
 * deck's HeatmapLayer kernel is `exp(-u²/0.0556) · 6/√π` with `u = r / radiusPixels`: a
 * Gaussian with σ = radiusPixels / 6 (screen pixels), unnormalised, so one point weighs
 * 6/√π at its centre. MDV's `contourFill` threshold and its contour spacing (a line every 2
 * weight units, see `smoothContour`) are in those weight units. DensityLayer's "absolute"
 * scale is weight per unit area, Σ w·exp(-r²/2σ²) / 2πσ², so
 *
 *     heatmapWeight = density · 2πσ² · 6/√π = density · 12√π · σ²
 *
 * which converts both exactly, rather than by eye.
 *
 * Deliberately not reproduced: the GLSL scales the fill's alpha, but not its colour, by
 * `min(value · intensity / (max · threshold), 1)`, so near zero density the fill adds to
 * the background rather than covering it. That glow is treated as a bug.
 *
 * MDV's radiusPixels is `30 · bandwidth · 2^zoom`, which is already zoom-constant in world
 * units: σ = 5 · bandwidth world units. DensityLayer's σ is in the data's own units, before
 * `modelMatrix`, so divide by the matrix's scale.
 */
const HEATMAP_WEIGHT_PER_DENSITY = 12 * Math.sqrt(Math.PI);
/** `contourFill` value MDV uses for "no fill" (see `contour_state.ts`). */
const NO_FILL = 10000;
/** One line every this many heatmap-weight units (MDV's `smoothContour`). */
const CONTOUR_SPACING = 2;
/**
 * The fill's edge strength relative to the fill, roughly as in MDV's GLSL contour shader
 * (`edgeShape · fillOpacity · 1.3`, alpha doubled). Tuned by eye, not for colour parity.
 * MDV has no separate edge setting yet; until it does, the edge follows the fill.
 */
const EDGE_PER_FILL = 2.6;

/** Geometric mean of the x and y basis lengths; 1 for a missing or degenerate matrix. */
function modelMatrixScale(matrix: Matrix4 | null | undefined): number {
    if (!matrix) return 1;
    const scale = Math.sqrt(Math.hypot(matrix[0], matrix[1], matrix[2]) * Math.hypot(matrix[4], matrix[5], matrix[6]));
    return Number.isFinite(scale) && scale > 0 ? scale : 1;
}

/**
 * Identifies what a `getWeight` computes, since the function itself is a new closure on
 * every render (field contours rebuild their props on each zoom step): an object compared
 * by identity (e.g. the field's data), then values compared by value (e.g. its min/max).
 */
export type WeightKey = readonly [object, ...unknown[]];

type WeightCacheEntry = { params: readonly unknown[]; length: number; weights: Float32Array };
const UNWEIGHTED = {};
/**
 * Keyed on the index list's identity, which is stable until the filter or category
 * changes, then on the weight key's object, so fields sharing one index list each keep
 * their own entry.
 */
const weightCache = new WeakMap<object, WeakMap<object, WeightCacheEntry>>();

/**
 * 1 (or the field weight) at the indices a contour draws from, 0 elsewhere.
 *
 * A field's `getWeight` is NaN where the row has no value (sparse gene scores normalise
 * `(value - min) / range`). The WebGL heatmap's blending happens to drop those; splatted
 * here, one NaN poisons every cell its kernel touches, so they are left out explicitly.
 */
function weightsFor(
    indices: ArrayLike<number>,
    length: number,
    getWeight?: (i: number) => number,
    weightKey?: WeightKey,
) {
    // without a key, a weighted call caches on the closure: right, but only while it lives
    const [keyObject, ...params]: WeightKey = getWeight ? (weightKey ?? [getWeight]) : [UNWEIGHTED];
    let byKey = weightCache.get(indices);
    if (!byKey) {
        byKey = new WeakMap();
        weightCache.set(indices, byKey);
    }
    const cached = byKey.get(keyObject);
    if (
        cached &&
        cached.length === length &&
        cached.params.length === params.length &&
        cached.params.every((p, k) => Object.is(p, params[k]))
    ) {
        return cached.weights;
    }
    const weights = new Float32Array(length);
    for (let k = 0; k < indices.length; k++) {
        const i = indices[k];
        const w = getWeight ? getWeight(i) : 1;
        weights[i] = Number.isFinite(w) ? w : 0;
    }
    byKey.set(keyObject, { params, length, weights });
    return weights;
}

export function contourToDensityProps(
    props: NonNullable<ContourLayerProps> & {
        getWeight?: (i: number) => number;
        weightKey?: WeightKey;
        modelMatrix?: Matrix4 | null;
    },
): DensityLayerProps {
    const { id, x, y, data, getWeight, weightKey, bandwidth, fillOpacity, contourOpacity, contourFill, colorRange } = props;
    // radiusPixels / 6 screen pixels, and `radiusPixels = 30 · bandwidth · 2^zoom`: 5 · bandwidth world units.
    const sigma = (5 * bandwidth) / modelMatrixScale(props.modelMatrix);
    const perDensity = HEATMAP_WEIGHT_PER_DENSITY * sigma * sigma;
    return {
        id,
        x,
        y,
        weights: weightsFor(data, x.length, getWeight, weightKey),
        bandwidth: sigma,
        scale: "absolute",
        colorRange: colorRange.map((c) => [c[0], c[1], c[2]] as const),
        fillLevel: contourFill >= NO_FILL ? null : contourFill / perDensity,
        fillOpacity,
        edgeOpacity: Math.min(1, EDGE_PER_FILL * fillOpacity),
        contourInterval: CONTOUR_SPACING / perDensity,
        contourOpacity,
        modelMatrix: props.modelMatrix ?? null,
        pickable: false,
    };
}
