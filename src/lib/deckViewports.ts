/** The parts of a deck viewport needed to find which one a canvas point is in. */
export type CanvasViewport = {
    id: string;
    x?: number;
    y?: number;
    width?: number;
    height?: number;
};

/**
 * The viewport containing the canvas point, for a canvas that draws several views (a
 * density grid). The last match wins, as it is drawn on top; off every viewport, the first.
 */
export function getViewportAtCanvasPoint<V extends CanvasViewport>(
    viewports: readonly V[],
    canvasX: number,
    canvasY: number,
): V | undefined {
    for (let i = viewports.length - 1; i >= 0; i--) {
        const { x = 0, y = 0, width = 0, height = 0 } = viewports[i];
        if (canvasX >= x && canvasX < x + width && canvasY >= y && canvasY < y + height) return viewports[i];
    }
    return viewports[0];
}
