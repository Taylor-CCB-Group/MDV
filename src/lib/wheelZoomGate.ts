import { OrbitController, OrthographicController } from "@deck.gl/core";
import type { MjolnirWheelEvent } from "mjolnir.js";

/**
 * Wheel zoom that shares the wheel with the page.
 *
 * When a chart sits in something that scrolls (the dashboard, or a grid of cells), a plain wheel
 * scrolls it and only Ctrl / ⌘ + wheel zooms the chart; the chart shows a hint saying so (see
 * `WheelZoomHint`). A trackpad pinch arrives as Ctrl + wheel, so it zooms either way. With nothing
 * to scroll (a chart popped out or full screen, or a dashboard that fits) a plain wheel zooms.
 */

/** Dispatched (bubbling) from the canvas when a plain wheel is left to scroll the page. */
export const WHEEL_ZOOM_PASSED_EVENT = "mdv-wheel-zoom-passed";

/** Wheel events closer together than this belong to one gesture, momentum included. */
const GESTURE_GAP_MS = 200;

const SCROLLING_OVERFLOW = new Set(["auto", "scroll", "overlay"]);

/** Whether something between `element` and the top of its document (or full-screen element) scrolls vertically. */
export function hasScrollableAncestor(element: Element | null): boolean {
    const doc = element?.ownerDocument;
    const view = doc?.defaultView;
    if (!doc || !view) return false;
    for (let el = element.parentElement; el; el = el.parentElement) {
        if (el.scrollHeight > el.clientHeight + 1) {
            const { overflowY } = view.getComputedStyle(el);
            // the root scrolls the viewport unless it hides its overflow
            if (el === doc.scrollingElement ? overflowY !== "hidden" : SCROLLING_OVERFLOW.has(overflowY)) {
                return true;
            }
        }
        if (el === doc.fullscreenElement) return false;
    }
    return false;
}

function hasZoomModifier(event: WheelEvent) {
    return event.ctrlKey || event.metaKey;
}

/**
 * Whether a wheel event zooms, decided once per gesture so a zoom keeps going when the modifier
 * is let go mid-swipe, and a scroll keeps scrolling when the pointer comes onto the chart.
 */
class WheelGestureGate {
    private lastTime = Number.NEGATIVE_INFINITY;
    private zooming = false;

    shouldZoom(event: MjolnirWheelEvent): boolean {
        const src = event.srcEvent;
        const now = src.timeStamp;
        const continuing = now - this.lastTime < GESTURE_GAP_MS;
        this.lastTime = now;
        // a modifier starts zooming at any point; without one, carry on as the gesture began
        if (hasZoomModifier(src)) this.zooming = true;
        else if (!continuing) {
            this.zooming = !hasScrollableAncestor(src.target instanceof Element ? src.target : null);
        }
        if (!this.zooming && !continuing && src.target instanceof EventTarget) {
            src.target.dispatchEvent(new CustomEvent(WHEEL_ZOOM_PASSED_EVENT, { bubbles: true }));
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
