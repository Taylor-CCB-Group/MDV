import { OrbitViewport } from "@deck.gl/core";

/**
 * The 3D scatter plot's camera: deck's `OrbitView`, and how the classic (WebGL) 3D scatter's
 * camera maps onto it.
 */

/** Which data axis points up on screen, and is the one the view orbits round. */
export type OrbitAxis = "Y" | "Z";

/** A data unit's length in world units along x, y and z (the classic 3D scatter's `axis_scales`). */
export type AxisScales = [number, number, number];

export type OrbitViewParams = {
    target: [number, number, number];
    zoom: number;
    rotationOrbit: number;
    rotationX: number;
};

/** The camera a 3D selection was drawn with, which says what its outline enclosed. */
export type OrbitSelectionView = OrbitViewParams & {
    width: number;
    height: number;
    orbitAxis: OrbitAxis;
};

/**
 * The view a fitted 3D chart starts from: a little above and to the side, which is where the
 * classic 3D scatter put its camera (`theta: -1.038, phi: 0.261`, see {@link orbitViewFromLegacyCamera}).
 */
export const DEFAULT_ROTATION_ORBIT = 30.5;
export const DEFAULT_ROTATION_X = 15;

/** Leaves room round the data when fitting a view to it. */
const FIT_PADDING_ZOOM = 0.3;

/**
 * The vertical field of view, in degrees: the classic 3D scatter's, so a converted chart has the
 * same perspective. (Zoom sets the scale at the target; this only sets how strong the perspective is.)
 */
export const ORBIT_FOVY = 45;
const LEGACY_FOVY = (ORBIT_FOVY * Math.PI) / 180;

export const DEFAULT_AXIS_SCALES: AxisScales = [1, 1, 1];

/** `axis_scales` if it is three positive numbers, else unscaled. */
export function validAxisScales(scales: unknown): AxisScales {
    if (Array.isArray(scales) && scales.length === 3 && scales.every((s) => typeof s === "number" && s > 0)) {
        return [scales[0], scales[1], scales[2]];
    }
    return DEFAULT_AXIS_SCALES;
}

/** The corners of a box, x, y and z. */
type Bounds3D = { min: readonly number[]; max: readonly number[] };

/**
 * The view that shows `bounds` (in world units, so already axis-scaled) whichever way it is turned:
 * the sphere round the box fills the shorter side of the chart.
 */
export function fitOrbitView(
    { min, max }: Bounds3D,
    width: number,
    height: number,
    rotation: Pick<OrbitViewParams, "rotationOrbit" | "rotationX"> = {
        rotationOrbit: DEFAULT_ROTATION_ORBIT,
        rotationX: DEFAULT_ROTATION_X,
    },
): OrbitViewParams {
    const mid = (i: number) => (min[i] + max[i]) / 2;
    const target: [number, number, number] = [mid(0), mid(1), mid(2)];
    const diameter = Math.hypot(...[0, 1, 2].map((i) => max[i] - min[i]));
    const side = Math.min(width, height);
    let zoom = Math.log2(side / Math.max(diameter, 1e-9)) - FIT_PADDING_ZOOM;
    if (!Number.isFinite(zoom)) zoom = 0;
    return { target, zoom, ...rotation };
}

type LegacyCamera = { distance: number; theta: number; phi: number };

/**
 * The deck view that shows what the classic 3D scatter's `camera` did, for a `orbitAxis: "Y"` view
 * whose points are scaled by `axis_scales`.
 *
 * The classic chart drew a point at `2 * (-x * sx, y * sy, -z * sz)` (a half turn about y of the
 * scaled data), and looked at `2 * (-cx, cy, -cz)` (`center`, which is already scaled) from
 * `distance` away, `theta` round from x towards z and `phi` up, with y up and a 45° field of view.
 * So the target is `center` itself; deck's `rotationX` is `phi`; its `rotationOrbit` is `theta` a
 * quarter turn on, deck starting from +z where the classic chart started from +x; and a data unit
 * at the target is `plotHeight / (distance * tan(fovy / 2))` pixels.
 * @param plotHeight the height of the classic chart's plot, in pixels
 */
export function orbitViewFromLegacyCamera(
    camera: LegacyCamera,
    center: [number, number, number],
    plotHeight: number,
): OrbitViewParams | undefined {
    const { distance, theta, phi } = camera;
    if (![distance, theta, phi, plotHeight, ...center].every(Number.isFinite) || distance <= 0 || plotHeight <= 0) {
        return undefined;
    }
    const degrees = (r: number) => (r * 180) / Math.PI;
    // deck keeps the orbit angle in [-180, 180)
    const orbit = ((((degrees(theta) + 90 + 180) % 360) + 360) % 360) - 180;
    return {
        target: [...center],
        zoom: Math.log2(plotHeight / (distance * Math.tan(LEGACY_FOVY / 2))),
        rotationOrbit: orbit,
        rotationX: degrees(phi),
    };
}

/** Within a pixel and a hundredth of a degree or zoom level: the view a selection was drawn in. */
export function isSameOrbitView(a: OrbitSelectionView, b: OrbitSelectionView) {
    const near = (x: number, y: number, tolerance: number) => Math.abs(x - y) <= tolerance;
    return (
        a.orbitAxis === b.orbitAxis &&
        near(a.width, b.width, 1) &&
        near(a.height, b.height, 1) &&
        near(a.zoom, b.zoom, 0.01) &&
        near(a.rotationOrbit, b.rotationOrbit, 0.01) &&
        near(a.rotationX, b.rotationX, 0.01) &&
        // a hundredth of a pixel's worth of data at this zoom
        a.target.every((t, i) => near(t, b.target[i], 0.01 / 2 ** a.zoom))
    );
}

/** The deck viewport for `view`, as the chart draws it. */
export function orbitViewport(view: OrbitSelectionView) {
    return new OrbitViewport({ ...view, x: 0, y: 0, fovy: ORBIT_FOVY });
}

/**
 * Where each row's point is drawn in `view`, in pixels from the top left of the plot; NaN for a
 * point behind the camera. `modelMatrix` takes data to world units (axis and region scaling).
 */
export function projectRowsToScreen(
    columns: [ArrayLike<number>, ArrayLike<number>, ArrayLike<number>],
    modelMatrix: ArrayLike<number>,
    view: OrbitSelectionView,
): [Float32Array, Float32Array] {
    const [xs, ys, zs] = columns;
    const n = xs.length;
    const px = new Float32Array(n);
    const py = new Float32Array(n);
    // the matrix deck projects with (`Viewport.project`), after the layer's model matrix; column-major
    const m = multiply(orbitViewport(view).pixelProjectionMatrix, modelMatrix);
    for (let i = 0; i < n; i++) {
        const x = xs[i];
        const y = ys[i];
        const z = zs[i];
        const w = m[3] * x + m[7] * y + m[11] * z + m[15];
        if (!(w > 0)) {
            px[i] = Number.NaN;
            py[i] = Number.NaN;
            continue;
        }
        px[i] = (m[0] * x + m[4] * y + m[8] * z + m[12]) / w;
        py[i] = (m[1] * x + m[5] * y + m[9] * z + m[13]) / w;
    }
    return [px, py];
}

/** `a * b`, both 4x4 column-major. */
function multiply(a: ArrayLike<number>, b: ArrayLike<number>) {
    const out = new Float64Array(16);
    for (let col = 0; col < 4; col++) {
        for (let row = 0; row < 4; row++) {
            let sum = 0;
            for (let k = 0; k < 4; k++) sum += a[k * 4 + row] * b[col * 4 + k];
            out[col * 4 + row] = sum;
        }
    }
    return out;
}

const isNumber = (x: unknown): x is number => typeof x === "number" && Number.isFinite(x);

/** The view saved with a 3D selection (in its feature's `properties.view`), if it has a valid one. */
export function orbitSelectionViewOf(properties: unknown): OrbitSelectionView | undefined {
    if (!properties || typeof properties !== "object" || !("view" in properties)) return undefined;
    const view: unknown = properties.view;
    if (!view || typeof view !== "object") return undefined;
    const v: Partial<Record<keyof OrbitSelectionView, unknown>> = view;
    const { target, zoom, rotationOrbit, rotationX, width, height, orbitAxis } = v;
    if (!Array.isArray(target) || target.length !== 3 || !target.every(isNumber)) return undefined;
    if (![zoom, rotationOrbit, rotationX, width, height].every(isNumber)) return undefined;
    if (orbitAxis !== "Y" && orbitAxis !== "Z") return undefined;
    return {
        target: [target[0], target[1], target[2]],
        zoom: Number(zoom),
        rotationOrbit: Number(rotationOrbit),
        rotationX: Number(rotationX),
        width: Number(width),
        height: Number(height),
        orbitAxis,
    };
}
