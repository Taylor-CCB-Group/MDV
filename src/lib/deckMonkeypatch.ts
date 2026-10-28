import type { EventManager } from 'mjolnir.js';
import { Deck } from '@deck.gl/core';
import { EditableGeoJsonLayer } from '@deck.gl-community/editable-layers';
import { cpuEditablePicks } from '@/editable-layers/cpuPicks';

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

export class MonkeyPatchEditableGeoJsonLayer extends EditableGeoJsonLayer {
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
     */
    getPicks(screenCoords: Parameters<EditableGeoJsonLayer["getPicks"]>[0]): EditablePicks {
        if (this.context.device.type !== "webgpu") return super.getPicks(screenCoords);
        const [x, y] = screenCoords;
        const guides = this.getActiveMode().getGuides(this.getModeProps(this.props));
        const picks = cpuEditablePicks({
            features: this.props.data.features,
            guides: guides?.features ?? [],
            point: screenCoords,
            project: (position) => this.project([position[0], position[1], 0]),
            radius: this.props.pickingRadius ?? 0,
        });
        return picks.map((pick) => ({ ...pick, color: null, layer: this, x, y, pixelRatio: 1 }));
    }
}

/** Deck's own (private) helpers, so the replacement manager is configured exactly as deck would. */
function getDeckPrivates(deck: MonkeyPatchDeck) {
    const createEventManager: unknown = Reflect.get(deck, '_createEventManager');
    const getEventRoot: unknown = Reflect.get(deck, '_getEventRoot');
    const isMultiCanvasMode: unknown = Reflect.get(deck, '_isMultiCanvasMode');
    if (typeof createEventManager !== 'function' || typeof getEventRoot !== 'function' || typeof isMultiCanvasMode !== 'function') {
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
    // experimental multi-canvas decks keep a manager per canvas; we don't use them.
    if (privates.isMultiCanvasMode()) return;

    const eventManager = privates.createEventManager(privates.getEventRoot(deck.canvas));
    deck.eventManager = eventManager;
    deck.eventManagers = Object.fromEntries(
        Object.entries(deck.eventManagers).map(([canvasId, em]) => [canvasId, em === oldEventManager ? eventManager : em]),
    );
    // the view manager prefers `eventManagers[canvasId]`, falling back to its own `_eventManager`;
    // with a new manager it finalizes the old controllers and builds new ones on the next update.
    Reflect.set(deck.viewManager, '_eventManager', eventManager);
    deck.viewManager.setProps({ eventManagers: deck.eventManagers });

    // Editable layers subscribe to the event manager once, in initializeState, so the ones already
    // in the deck are still attached to the old one. Move their existing handlers across; layers
    // initialized after this subscribe to the new manager themselves.
    // These are looked up in the deck rather than passed in: the instance a chart created may not
    // be the one in the deck (SpatialData.js re-ids host layers by cloning them).
    const layers: unknown[] = deck.layerManager?.getLayers() ?? [];
    for (const layer of layers) {
        if (layer instanceof EditableGeoJsonLayer) layer._addEventHandlers();
    }

    // only now - controllers unsubscribe from the old manager as they are finalized above.
    oldEventManager.destroy();
}
