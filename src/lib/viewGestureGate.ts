import { OrbitController, OrthographicController } from "@deck.gl/core";
import type { MjolnirGestureEvent, MjolnirWheelEvent } from "mjolnir.js";

/**
 * Which pointer gestures move a chart's view, for the deck charts with a selection toolbar.
 *
 * - **Wheel**: Ctrl / ⌘ + wheel zooms the chart, and a plain wheel is left to scroll whatever is
 *   around it. A trackpad pinch arrives as Ctrl + wheel, so it zooms too.
 * - **Drag**: a plain drag draws a selection with the current tool; ⌘-drag (Ctrl-drag off a Mac)
 *   or a right- or middle-button drag moves the view: pans in 2D, orbits in 3D, and pans in 3D
 *   with Shift as well. The Pan tool is a lock that makes a plain drag move the view too.
 * - **Touch**: one finger draws; two fingers pan (orbit in 3D) and pinch to zoom.
 *
 * This holds wherever the chart is (dashboard, grid of cells, popout, full screen), so there is
 * one habit to learn. A hint says how to zoom or move when a plain gesture doesn't (see
 * `ViewGestureHint`).
 */

/** Dispatched (bubbling) from the canvas when a plain wheel is left to scroll the page. */
export const WHEEL_ZOOM_PASSED_EVENT = "mdv-wheel-zoom-passed";

/** Dispatched (bubbling) from the canvas when a plain drag draws a selection rather than moving the view. */
export const VIEW_DRAG_PASSED_EVENT = "mdv-view-drag-passed";

export const isMac = typeof navigator !== "undefined" && /Mac|iPhone|iPad/.test(navigator.userAgent);

/** ⌘ on a Mac, where Ctrl-click is a right-click; Ctrl elsewhere. */
export const VIEW_MODIFIER_NAME = isMac ? "⌘" : "Ctrl";

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

type ModifierState = { metaKey?: boolean; ctrlKey?: boolean; shiftKey?: boolean };

/** Whether the key that moves the view is held (see {@link VIEW_MODIFIER_NAME}). */
export function hasViewModifier(src: ModifierState, mac = isMac) {
    return mac ? Boolean(src.metaKey) : Boolean(src.ctrlKey || src.metaKey);
}

/** The modifier keys and buttons of the event a gesture started with, where it has them. */
function sourceOf(event: Pick<MjolnirGestureEvent, "srcEvent">): ModifierState & { buttons?: number } {
    const src: unknown = event.srcEvent;
    return src && typeof src === "object" ? src : {};
}

/**
 * Whether a drag moves the view whatever tool is chosen: with the view modifier held, or with the
 * right or middle button.
 */
export function isViewDrag(event: Pick<MjolnirGestureEvent, "srcEvent" | "rightButton" | "middleButton">, mac = isMac) {
    // `buttons` of the move that started the drag: 2 is the right button, 4 the middle one
    const buttons = sourceOf(event).buttons ?? 0;
    return Boolean(event.rightButton || event.middleButton || buttons & 6 || hasViewModifier(sourceOf(event), mac));
}

const LEARNED_KEY = "mdv.viewDragLearned";

/** Whether this browser has moved a view with a drag that needed the modifier (or another button). */
export function hasLearnedViewDrag() {
    try {
        return localStorage.getItem(LEARNED_KEY) === "1";
    } catch {
        return false;
    }
}

function learnViewDrag() {
    try {
        localStorage.setItem(LEARNED_KEY, "1");
    } catch {
        // without storage the hint keeps showing, which is harmless
    }
}

/** Props the gated controllers take on top of deck's own (`controller: {type, ...}`). */
export type ViewDragProps = {
    /** The Pan tool's lock: a plain drag moves the view too. */
    dragMovesView?: boolean;
};

/**
 * How a drag that starts now moves the view, or null if it leaves the view alone (and draws).
 * `orbit`: the view turns about its target, and Shift pans it instead.
 */
function viewDragMode(event: MjolnirGestureEvent, dragMovesView: boolean, orbit: boolean) {
    const explicit = isViewDrag(event);
    if (!explicit && !dragMovesView) return null;
    if (explicit) learnViewDrag();
    return orbit && !sourceOf(event).shiftKey ? "rotate" : "pan";
}

/** The event as a left-button drag, so deck doesn't read the button as a choice of mode. */
function asLeftButton(event: MjolnirGestureEvent): MjolnirGestureEvent {
    return { ...event, rightButton: false, middleButton: false };
}

/*
 * Deck chooses between its two drag modes by whether any modifier key or the right button is
 * held (pan for 3D, rotate for 2D, which a 2D view can't), inverted when `dragMode` is "pan".
 * The gated controllers choose the mode themselves (`viewDragMode`), and answer deck's question
 * about modifier keys so that it starts that mode.
 */

export class GatedOrthographicController extends OrthographicController {
    private wheelGate = new WheelGestureGate();
    private dragMovesView = false;
    private functionKeyOverride: boolean | null = null;

    setProps(props: Parameters<OrthographicController["setProps"]>[0] & ViewDragProps) {
        this.dragMovesView = Boolean(props.dragMovesView);
        super.setProps(props);
    }

    isFunctionKeyPressed(event: Parameters<OrthographicController["isFunctionKeyPressed"]>[0]) {
        return this.functionKeyOverride ?? super.isFunctionKeyPressed(event);
    }

    protected _onWheel(event: MjolnirWheelEvent): boolean {
        if (!this.isPointInBounds(this.getCenter(event), event)) return false;
        return this.wheelGate.shouldZoom(event) && super._onWheel(event);
    }

    protected _onPanStart(event: MjolnirGestureEvent): boolean {
        const mode = viewDragMode(event, this.dragMovesView, false);
        if (!mode) return false;
        this.functionKeyOverride = (mode === "pan") !== (this.dragMode === "pan");
        try {
            return super._onPanStart(asLeftButton(event));
        } finally {
            this.functionKeyOverride = null;
        }
    }
}

export class GatedOrbitController extends OrbitController {
    private wheelGate = new WheelGestureGate();
    private dragMovesView = false;
    private functionKeyOverride: boolean | null = null;

    setProps(props: Parameters<OrbitController["setProps"]>[0] & ViewDragProps) {
        this.dragMovesView = Boolean(props.dragMovesView);
        super.setProps(props);
    }

    isFunctionKeyPressed(event: Parameters<OrbitController["isFunctionKeyPressed"]>[0]) {
        return this.functionKeyOverride ?? super.isFunctionKeyPressed(event);
    }

    protected _onWheel(event: MjolnirWheelEvent): boolean {
        if (!this.isPointInBounds(this.getCenter(event), event)) return false;
        return this.wheelGate.shouldZoom(event) && super._onWheel(event);
    }

    protected _onPanStart(event: MjolnirGestureEvent): boolean {
        const mode = viewDragMode(event, this.dragMovesView, true);
        if (!mode) return false;
        this.functionKeyOverride = (mode === "pan") !== (this.dragMode === "pan");
        try {
            return super._onPanStart(asLeftButton(event));
        } finally {
            this.functionKeyOverride = null;
        }
    }
}
