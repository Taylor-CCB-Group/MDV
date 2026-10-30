/**
 * How a 2D scatter plot maps data units to pixels.
 * - `equal`: one zoom for both axes, so a unit is the same length in x and y. Right for physical
 *   coordinates, and what charts saved before this setting existed do.
 * - `fit`: each axis zooms separately, so data with very different x and y ranges fills the chart
 *   rather than collapsing along the axis with the smaller range. Points stay circular.
 */
export type AxisScaling = "equal" | "fit";

export type DataBounds = { minX: number; maxX: number; minY: number; maxY: number };

/** Leaves room round the data when fitting a view to it. */
const FIT_PADDING_ZOOM = 0.6;

export type OrthographicFit = {
    target: [number, number, number];
    /** The smaller of the axis zooms, which is what deck sizes points by. */
    zoom: number;
    /** Present only for `fit`. */
    zoomX?: number;
    zoomY?: number;
};

/** The view that shows `bounds` in a `width` x `height` chart. */
export function fitOrthographicView(
    bounds: DataBounds,
    width: number,
    height: number,
    scaling: AxisScaling,
): OrthographicFit {
    const { minX, maxX, minY, maxY } = bounds;
    // a relative epsilon, so that a single point or a constant column still gives a finite zoom
    const epsilon = Math.max(1e-9, Math.min(Math.abs(minX), Math.abs(maxX)) * 1e-6);
    const dx = Math.max(epsilon, maxX - minX);
    const dy = Math.max(epsilon, maxY - minY);
    const target: [number, number, number] = [(minX + maxX) / 2, (minY + maxY) / 2, 0];
    const zoomX = Math.log2(width / dx) - FIT_PADDING_ZOOM;
    const zoomY = Math.log2(height / dy) - FIT_PADDING_ZOOM;
    const zoom = Math.min(zoomX, zoomY);
    if (!Number.isFinite(zoom)) {
        console.warn("Invalid zoom value detected, using fallback");
        return { target, zoom: 0 };
    }
    return scaling === "fit" ? { target, zoom, zoomX, zoomY } : { target, zoom };
}

/**
 * The data range that point radius is proportional to. With `equal` that is the larger axis range.
 * With `fit` the axes are drawn at different scales, so it is the range of a square of the chart's
 * shorter side, which makes points about the same size on screen however the data is distributed.
 */
export function radiusReferenceRange(
    xRange: number,
    yRange: number,
    chartSize: { width: number; height: number },
    scaling: AxisScaling,
) {
    const { width, height } = chartSize;
    if (scaling !== "fit" || !(width > 0) || !(height > 0)) return Math.max(xRange, yRange);
    return Math.max(xRange / width, yRange / height) * Math.min(width, height);
}

type ZoomState = { zoom?: number | [number, number]; zoomX?: number; zoomY?: number };

/** A view state's zoom reduced to the plain numbers we save (deck reports `zoom` as a pair when the axes differ). */
export function serialisableZoom(viewState: ZoomState, scaling: AxisScaling) {
    const { zoom, zoomX, zoomY } = viewState;
    const scalar = Array.isArray(zoom) ? Math.min(...zoom) : (zoom ?? 0);
    if (scaling !== "fit" || zoomX === undefined || zoomY === undefined) return { zoom: scalar };
    return { zoom: Math.min(scalar, zoomX, zoomY), zoomX, zoomY };
}

/**
 * The classic scatter's `radius` is a point's diameter in pixels at the initial view; the deck chart's
 * `radius * course_radius` is a fraction (1/10000) of the data range, so the pixels depend on the chart.
 * This is the deck setting that draws the same points in a chart of `size`, assuming the data is about
 * as wide as it is tall (which is not known when a config is converted). The slider tops out at 20, so
 * larger values use `course_radius` instead.
 */
export function deckRadiusForClassicRadius(classicRadius: number, size: [number, number]) {
    const plotSide = Math.min(...size);
    // at the fitted view a deck radius of r is about r / 10000 * 2^-FIT_PADDING_ZOOM * plotSide pixels
    const total = ((classicRadius / 2) * 10000) / (2 ** -FIT_PADDING_ZOOM * plotSide);
    if (!(total > 0) || !Number.isFinite(total)) return undefined;
    if (total <= 20) return { radius: Math.round(total * 10) / 10, course_radius: 1 };
    return { radius: Math.round(total) / 10, course_radius: 10 };
}
