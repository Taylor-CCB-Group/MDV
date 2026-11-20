import { OrbitController, OrthographicController } from "@deck.gl/core";
import type { MjolnirWheelEvent } from "mjolnir.js";

/**
 * Wheel zoom behind a modifier: Ctrl / ⌘ + wheel zooms the chart, and a plain wheel is left to
 * scroll whatever is around it, with a hint saying how to zoom (see `WheelZoomHint`). A trackpad
 * pinch arrives as Ctrl + wheel, so it zooms too. This holds wherever the chart is (dashboard,
 * grid of cells, popout, full screen), so there is one habit to learn.
 */

/** Dispatched (bubbling) from the canvas when a plain wheel is left to scroll the page. */
export const WHEEL_ZOOM_PASSED_EVENT = "mdv-wheel-zoom-passed";

/** Wheel events closer together than this belong to one gesture, momentum included. */
const GESTURE_GAP_MS = 200;

/**
 * Whether a wheel event zooms, decided once per gesture so a zoom keeps going when the modifier
 * is let go mid-swipe, rather than the rest of the swipe scrolling the page.
 */
class WheelGestureGate {
    private lastTime = Number.NEGATIVE_INFINITY;
    private zooming = false;

    shouldZoom(event: MjolnirWheelEvent): boolean {
        const src = event.srcEvent;
        // this window's clock: event time stamps count from their own window's start, which
        // differs in a popout
        const now = performance.now();
        const continuing = now - this.lastTime < GESTURE_GAP_MS;
        this.lastTime = now;
        if (src.ctrlKey || src.metaKey) this.zooming = true;
        else if (!continuing) {
            this.zooming = false;
            src.target?.dispatchEvent(new CustomEvent(WHEEL_ZOOM_PASSED_EVENT, { bubbles: true }));
        }
        return this.zooming;
    }
}

export class GatedOrthographicController extends OrthographicController {
    private wheelGate = new WheelGestureGate();

    protected _onWheel(event: MjolnirWheelEvent): boolean {
        if (!this.isPointInBounds(this.getCenter(event), event)) return false;
        return this.wheelGate.shouldZoom(event) && super._onWheel(event);
    }
}

export class GatedOrbitController extends OrbitController {
    private wheelGate = new WheelGestureGate();

    protected _onWheel(event: MjolnirWheelEvent): boolean {
        if (!this.isPointInBounds(this.getCenter(event), event)) return false;
        return this.wheelGate.shouldZoom(event) && super._onWheel(event);
    }
}
