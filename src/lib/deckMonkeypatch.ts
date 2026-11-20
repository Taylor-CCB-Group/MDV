import { Deck, type Viewport } from "@deck.gl/core";
import { EditableGeoJsonLayer } from "@deck.gl-community/editable-layers";
import type { EventManager } from "mjolnir.js";
import { cpuEditablePicks } from "@/editable-layers/cpuPicks";
import { getViewportAtCanvasPoint } from "./deckViewports";

// we aren't really using this as an actual subclass, just declaring things as public so we can access them.
//! we don't want to keep this around for future maintenance, hopefully deck can be patched so this is less necessary
// or we find a better way to fix things.
class MonkeyPatchDeck extends Deck<any> {
    declare public eventManager: EventManager | null;
    declare public eventManagers: Record<string, EventManager>;
    declare public viewManager;
    declare public layerManager;
    declare public canvas: HTMLCanvasElement | null;
}

type EditablePicks = ReturnType<EditableGeoJsonLayer["getPicks"]>;
type MjolnirEvent = Parameters<EditableGeoJsonLayer["toBasePointerEvent"]>[0];

/** The mjolnir events editable-layers subscribes to (its own `EVENT_TYPES`, which it does not export). */
const EDITABLE_EVENT_TYPES = ["click", "pointermove", "panstart", "panmove", "panend", "keyup", "dblclick"];

/** The same deck, with the members we reach for declared public. */
function asPatchDeck(deck: Deck<any>): MonkeyPatchDeck {
    return deck as MonkeyPatchDeck;
}

/**
 * Every event manager a deck listens with: one per canvas in deck's experimental multi-canvas
 * mode (`_canvases`, as the density grid uses), otherwise just the one.
 */
function deckEventManagers(deck: MonkeyPatchDeck): EventManager[] {
    const managers = Object.values(deck.eventManagers ?? {});
    if (managers.length > 0) return managers;
    return deck.eventManager ? [deck.eventManager] : [];
}

/** Which of a multi-canvas deck's canvases an event came from (deck's own private lookup). */
function canvasIdOfEvent(deck: MonkeyPatchDeck, event: MjolnirEvent): string | undefined {
    const lookup: unknown = Reflect.get(deck, "_getCanvasIdFromEvent");
    if (typeof lookup !== "function") return undefined;
    const id: unknown = lookup.call(deck, event);
    return typeof id === "string" ? id : undefined;
}

/**
 * The managers each editable layer's handler is on (editable-layers itself only uses deck's
 * first). Keyed by the handler, which lives in layer state: deck makes a new layer instance
 * every render and hands it the state, so a field on the instance would forget.
 */
const subscribedManagers = new WeakMap<object, EventManager[]>();

export class MonkeyPatchEditableGeoJsonLayer extends EditableGeoJsonLayer {
    /** The canvas the event being handled came from, in a multi-canvas deck. */
    private eventCanvasId: string | undefined;

    private get patchDeck(): MonkeyPatchDeck | undefined {
        const deck = this.context.deck;
        return deck ? asPatchDeck(deck) : undefined;
    }

    /**
     * Listen on every canvas, so a selection can be drawn in any cell of a density grid. deck
     * adds and drops canvases as cells scroll into view: `syncEventHandlers` follows that.
     */
    _addEventHandlers() {
        const deck = this.patchDeck;
        if (!deck) return;
        const { eventHandler } = this.state._editableLayerState;
        const subscribed = subscribedManagers.get(eventHandler) ?? [];
        for (const manager of deckEventManagers(deck)) {
            if (subscribed.includes(manager)) continue;
            for (const type of EDITABLE_EVENT_TYPES) {
                // a higher priority than deck's controllers, so a drag can stop the view panning
                manager.on(type, eventHandler, { priority: 100 });
            }
            subscribed.push(manager);
        }
        subscribedManagers.set(eventHandler, subscribed);
    }

    _removeEventHandlers() {
        const { eventHandler } = this.state._editableLayerState;
        for (const manager of subscribedManagers.get(eventHandler) ?? []) {
            for (const type of EDITABLE_EVENT_TYPES) manager.off(type, eventHandler);
        }
        subscribedManagers.delete(eventHandler);
    }

    /** Follow deck's current event managers, after its canvases change. */
    syncEventHandlers() {
        const deck = this.patchDeck;
        if (!deck) return;
        const current = deckEventManagers(deck);
        const subscribed = subscribedManagers.get(this.state._editableLayerState.eventHandler) ?? [];
        if (current.length === subscribed.length && current.every((m) => subscribed.includes(m))) return;
        // deck destroys the managers of canvases it drops; unsubscribing from those is harmless
        this._removeEventHandlers();
        this._addEventHandlers();
    }

    toBasePointerEvent(event: MjolnirEvent) {
        const deck = this.patchDeck;
        this.eventCanvasId = deck ? canvasIdOfEvent(deck, event) : undefined;
        return super.toBasePointerEvent(event);
    }

    /**
     * The view under a point in the canvas the current event came from. editable-layers
     * projects through `context.viewport`, which is whichever view deck drew last, so a deck
     * with several views (one per canvas, or several in one) must find the right one itself.
     * Null for a deck with a single view, which the default handling suits.
     */
    private viewportAt([x, y]: readonly number[]): Viewport | null {
        const deck = this.patchDeck;
        if (!deck) return null;
        const canvasId = this.eventCanvasId;
        // the views of that canvas under the point; the last is drawn on top
        if (canvasId) return deck.getViewports({ x, y, canvasId }).at(-1) ?? null;
        const viewports = deck.getViewports();
        return viewports.length > 1 ? (getViewportAtCanvasPoint(viewports, x, y) ?? null) : null;
    }

    getMapCoords(screenCoords: Parameters<EditableGeoJsonLayer["getMapCoords"]>[0]) {
        const viewport = this.viewportAt(screenCoords);
        if (!viewport) return super.getMapCoords(screenCoords);
        return viewport.unproject([screenCoords[0] - viewport.x, screenCoords[1] - viewport.y]);
    }

    /**
     * editable-layers picks synchronously on every pointer event (`pickMultipleObjects`),
     * which deck.gl cannot do on a WebGPU device: it throws, and the throw aborted the
     * click before deck dispatched it to the other layers.
     *
     * On WebGPU, compute the picks on the CPU instead (`cpuEditablePicks`): the layer holds
     * one shape and its handles, so this is cheap and exact, and answers for where the
     * pointer is now. An async GPU pick could only answer for where it was, so a quick move
     * then press drew a new shape rather than moving the existing one. It also sidesteps two
     * deck 9.4 WebGPU picking bugs met along the way: rectangle picks throw under
     * @spatialdata's picking fix, and the outline's path sublayer returns garbage indices.
     * Remove once editable-layers picks asynchronously.
     *
     * On WebGL, pick in the canvas the event came from, in a multi-canvas deck.
     */
    getPicks(screenCoords: Parameters<EditableGeoJsonLayer["getPicks"]>[0]): EditablePicks {
        const [x, y] = screenCoords;
        if (this.context.device.type !== "webgpu") {
            const deck = this.patchDeck;
            const canvasId = this.eventCanvasId;
            if (!deck || !canvasId) return super.getPicks(screenCoords);
            return deck.pickMultipleObjects({
                x,
                y,
                canvasId,
                layerIds: [this.props.id],
                radius: this.props.pickingRadius,
                depth: this.props.pickingDepth,
            });
        }
        const guides = this.getActiveMode().getGuides(this.getModeProps(this.props));
        const picks = cpuEditablePicks({
            features: this.props.data.features,
            guides: guides?.features ?? [],
            point: screenCoords,
            project: this.projectAt(screenCoords),
            radius: this.props.pickingRadius ?? 0,
        });
        return picks.map((pick) => ({ ...pick, color: null, layer: this, x, y, pixelRatio: 1 }));
    }

    /** Projects to canvas pixels through the view under the pointer (see `viewportAt`). */
    private projectAt(screenCoords: readonly number[]) {
        const viewport = this.viewportAt(screenCoords);
        if (!viewport) return (position: readonly number[]) => this.project([position[0], position[1], 0]);
        return (position: readonly number[]) => {
            const [px, py, pz] = viewport.project([position[0], position[1], 0]);
            return [px + viewport.x, py + viewport.y, pz];
        };
    }
}

/** Point a deck's editable layers at its current event managers, after its canvases change. */
export function syncEditableLayerEvents(deck: Deck<any>) {
    const layers: unknown[] = asPatchDeck(deck).layerManager?.getLayers() ?? [];
    for (const layer of layers) {
        if (layer instanceof MonkeyPatchEditableGeoJsonLayer) layer.syncEventHandlers();
    }
}

/** Deck's own (private) helpers, so the replacement manager is configured exactly as deck would. */
function getDeckPrivates(deck: MonkeyPatchDeck) {
    const createEventManager: unknown = Reflect.get(deck, "_createEventManager");
    const getEventRoot: unknown = Reflect.get(deck, "_getEventRoot");
    const isMultiCanvasMode: unknown = Reflect.get(deck, "_isMultiCanvasMode");
    if (
        typeof createEventManager !== "function" ||
        typeof getEventRoot !== "function" ||
        typeof isMultiCanvasMode !== "function"
    ) {
        return null;
    }
    return {
        createEventManager: (root: HTMLElement): EventManager => createEventManager.call(deck, root),
        getEventRoot: (canvas: HTMLCanvasElement): HTMLElement => getEventRoot.call(deck, canvas),
        isMultiCanvasMode: (): boolean => isMultiCanvasMode.call(deck),
    };
}

/**
 * This is undesirable and liable to break in the future, but we need to monkeypatch
 * deck.gl to fix issues with event handling in popout windows.
 *
 * mjolnir binds pointermove / pointerup to the window that owned the event root when the
 * EventManager was created, and an EventManager's element can't be reassigned
 * (https://visgl.github.io/mjolnir.js/docs/api-reference/event-manager), so when a chart moves
 * to another window we replace the manager everywhere deck holds it.
 *
 * Apply this function to a deck instance whenever the outer container changes.
 * Written against deck.gl 9.4 - check it whenever deck is updated.
 */
export function rebindMouseEvents(deckO: Deck<any>) {
    // the one place we reach past deck's protected API
    const deck = deckO as MonkeyPatchDeck;
    const oldEventManager = deck.eventManager;
    if (!oldEventManager || !deck.canvas || !deck.viewManager) {
        // deck initializes its device asynchronously; it will create its manager in the right window.
        return;
    }
    const privates = getDeckPrivates(deck);
    if (!privates) {
        console.error("deckMonkeypatch: deck internals changed - popped-out charts may not receive mouse events");
        return;
    }
    if (privates.isMultiCanvasMode()) {
        rebindCanvasEventManagers(deck, privates);
        return;
    }

    const eventManager = privates.createEventManager(privates.getEventRoot(deck.canvas));
    deck.eventManager = eventManager;
    deck.eventManagers = Object.fromEntries(
        Object.entries(deck.eventManagers).map(([canvasId, em]) => [
            canvasId,
            em === oldEventManager ? eventManager : em,
        ]),
    );
    // the view manager prefers `eventManagers[canvasId]`, falling back to its own `_eventManager`;
    // with a new manager it finalizes the old controllers and builds new ones on the next update.
    Reflect.set(deck.viewManager, "_eventManager", eventManager);
    deck.viewManager.setProps({ eventManagers: deck.eventManagers });

    // Editable layers subscribe to the event manager once, in initializeState, so the ones already
    // in the deck are still attached to the old one. Move their existing handlers across; layers
    // initialized after this subscribe to the new manager themselves.
    // These are looked up in the deck rather than passed in: the instance a chart created may not
    // be the one in the deck (SpatialData.js re-ids host layers by cloning them).
    moveEditableLayerHandlers(deck);

    // only now - controllers unsubscribe from the old manager as they are finalized above.
    oldEventManager.destroy();
}

/**
 * Editable layers subscribe to the event manager once, in initializeState, so the ones already
 * in the deck are still attached to the old one. Move their existing handlers across; layers
 * initialized after this subscribe to the new manager themselves.
 * These are looked up in the deck rather than passed in: the instance a chart created may not
 * be the one in the deck (SpatialData.js re-ids host layers by cloning them).
 */
function moveEditableLayerHandlers(deck: MonkeyPatchDeck) {
    const layers: unknown[] = deck.layerManager?.getLayers() ?? [];
    for (const layer of layers) {
        if (layer instanceof MonkeyPatchEditableGeoJsonLayer) layer.syncEventHandlers();
        else if (layer instanceof EditableGeoJsonLayer) layer._addEventHandlers();
    }
}

/** The parts of deck's private `CanvasManager` (multi-canvas mode) that rebinding touches. */
type CanvasTargets = {
    targets: Record<string, { eventRoot: HTMLElement; eventManager: EventManager }>;
    order: string[];
    eventManagers: Record<string, EventManager>;
};

function isCanvasTargets(value: unknown): value is CanvasTargets {
    if (typeof value !== "object" || value === null) return false;
    const targets: unknown = Reflect.get(value, "targets");
    return (
        typeof targets === "object" &&
        targets !== null &&
        Array.isArray(Reflect.get(value, "order")) &&
        Object.values(targets).every(
            (t: unknown) =>
                typeof t === "object" &&
                t !== null &&
                Reflect.get(t, "eventRoot") instanceof Element &&
                typeof Reflect.get(t, "eventManager") === "object",
        )
    );
}

/**
 * `rebindMouseEvents` for a deck in experimental multi-canvas mode (`_canvases`), which keeps an
 * event manager per canvas in its private `CanvasManager`: replace each one, as for one canvas.
 * Canvases deck adds later get managers made in the window they are in by then.
 */
function rebindCanvasEventManagers(deck: MonkeyPatchDeck, privates: NonNullable<ReturnType<typeof getDeckPrivates>>) {
    const canvasManager: unknown = Reflect.get(deck, "_canvasManager");
    const viewManager = deck.viewManager;
    if (!isCanvasTargets(canvasManager) || !viewManager) {
        console.error("deckMonkeypatch: deck internals changed - popped-out grids may not receive mouse events");
        return;
    }
    const oldManagers: EventManager[] = [];
    for (const target of Object.values(canvasManager.targets)) {
        oldManagers.push(target.eventManager);
        target.eventManager = privates.createEventManager(target.eventRoot);
    }
    canvasManager.eventManagers = Object.fromEntries(
        Object.entries(canvasManager.targets).map(([id, target]) => [id, target.eventManager]),
    );
    deck.eventManagers = canvasManager.eventManagers;
    deck.eventManager = canvasManager.eventManagers[canvasManager.order[0]] ?? null;
    Reflect.set(viewManager, "_eventManager", deck.eventManager);
    viewManager.setProps({ eventManagers: deck.eventManagers });
    moveEditableLayerHandlers(deck);
    for (const manager of oldManagers) manager.destroy();
}
