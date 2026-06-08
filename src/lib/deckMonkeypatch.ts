import { EventManager, InputDirection, Pan, Pinch, Tap } from 'mjolnir.js';
// import { EVENT_HANDLERS, RECOGNIZERS } from '@deck.gl/core/dist/lib/constants';
import { Deck } from '@deck.gl/core';
import { EditableGeoJsonLayer } from '@deck.gl-community/editable-layers';
import type { Position } from '@turf/helpers';
import { unprojectCanvasPoint } from '@/react/components/densityGridUtils';

type DeckLayerLike = {
    id: string;
    context?: { deck?: unknown };
    getLayers?: () => unknown;
};

function flattenDeckLayers(layers: readonly DeckLayerLike[]): DeckLayerLike[] {
    const flattened: DeckLayerLike[] = [];
    for (const layer of layers) {
        flattened.push(layer);
        const sublayers = layer.getLayers?.();
        if (Array.isArray(sublayers)) {
            flattened.push(...flattenDeckLayers(sublayers as DeckLayerLike[]));
        }
    }
    return flattened;
}

/** Resolve the selection layer instance that is actually mounted on the deck. */
export function resolveDeckSelectionLayer(
    deck: Deck,
    selectionLayer?: EditableGeoJsonLayer,
): EditableGeoJsonLayer | undefined {
    const layerManager = (deck as unknown as { layerManager?: { layers?: DeckLayerLike[] } }).layerManager;
    const mounted = flattenDeckLayers(layerManager?.layers ?? []).find((layer) =>
        layer.id.startsWith("selection_"),
    );
    if (mounted) {
        return mounted as EditableGeoJsonLayer;
    }
    if (selectionLayer?.context?.deck) {
        return selectionLayer;
    }
    return undefined;
}

// we don't need the keys; we'll break in.
// we aren't really using this as an actual subclass, just declaring things as public so we can access them.
//! we don't want to keep this around for future maintenance, hopefully deck can be patched so this is less necessary
// or we find a better way to fix things.
class MonkeyPatchDeck extends Deck<any> {
    declare public eventManager;
    declare public viewManager;
    declare public animationLoop;
    declare public canvas;
    declare _onPointerDown: (event: unknown) => void;
    declare _onPointerMove: (event: unknown) => void;
    declare _onPointerUp?: (event: unknown) => void;
}

/**
In editable-layers.js, click and dblclick events weren't being handled possibly because of a dependency mismatch?

This should be removed when we next update deck.gl et al.

```typescript
const EVENT_TYPES = ['anyclick', 'pointermove', 'panstart', 'panmove', 'panend', 'keyup'];

class EditableLayer {
    _addEventHandlers() {
        // @ts-expect-error accessing protected props
        const { eventManager } = this.context.deck;
        const { eventHandler } = this.state._editableLayerState;
        for (const eventType of EVENT_TYPES) {
            eventManager.on(eventType, eventHandler, {
                // give nebula a higher priority so that it can stop propagation to deck.gl's map panning handlers
                priority: 100
            });
        }
    }
    _forwardEventToCurrentLayer(event) {
        const currentLayer = this.getCurrentLayer();
        // Use a naming convention to find the event handling function for this event type
        const func = currentLayer[`_on${event.type}`].bind(currentLayer);
        //...
    }
}
```
*/

/**
 * Editable layers use `context.viewport.unproject` with canvas-relative pixels.
 * In density-grid mode one canvas hosts many sub-viewports, so we resolve the
 * viewport under the pointer before unprojecting.
 */
export class MonkeyPatchEditableGeoJsonLayer extends EditableGeoJsonLayer {
    getMapCoords(screenCoords: Position): Position {
        const deck = this.context.deck;
        const viewports = deck?.getViewports?.();
        if (viewports && viewports.length > 1) {
            return unprojectCanvasPoint(viewports, screenCoords[0], screenCoords[1]) as Position;
        }
        return this.context.viewport.unproject([screenCoords[0], screenCoords[1]]) as Position;
    }

    getPicks(screenCoords: [number, number]) {
        const deck = this.context.deck;
        if (!deck) {
            return [];
        }
        const radius = this.props.pickingRadius;
        const depth = this.props.pickingDepth;
        const layerId = this.props.id;
        const allPicks = deck.pickMultipleObjects({
            x: screenCoords[0],
            y: screenCoords[1],
            radius,
            depth,
        });
        const matching = allPicks.filter(
            (pick) =>
                pick.layer?.id &&
                (pick.layer.id === layerId || pick.layer.id.startsWith(`${layerId}-`)),
        );
        if (matching.length > 0) {
            return matching;
        }
        return deck.pickMultipleObjects({
            x: screenCoords[0],
            y: screenCoords[1],
            layerIds: [layerId],
            radius,
            depth,
        });
    }
}

const EVENT_HANDLERS: { [eventName: string]: string } = {
    click: 'onClick', //is onClick a thing?
    anyclick: 'onClick',
    panstart: 'onDragStart',
    panmove: 'onDrag',
    panend: 'onDragEnd'
} as const;

const RECOGNIZERS = {
    multipan: [Pan, { threshold: 10, direction: InputDirection.Vertical, pointers: 2 }],
    pinch: [Pinch, {}, null, ['multipan']],
    pan: [Pan, { threshold: 1 }, ['pinch'], ['multipan']],
    dblclick: [Tap, { event: 'dblclick', taps: 2 }],
    // pointerdown is being handled by deck.js _onPointerDown which does picking...
    // and then we miss out.
    click: [Tap, { event: 'click' }, null, ['dblclick']],
} as const;

const DECK_EVENT_MANAGER_EVENTS = {
    pointerdown: (deck: MonkeyPatchDeck) => deck._onPointerDown,
    pointermove: (deck: MonkeyPatchDeck) => deck._onPointerMove,
    pointerleave: (deck: MonkeyPatchDeck) => deck._onPointerMove,
    // Without pointerup, tap/click recognizers and panend do not reliably fire.
    pointerup: (deck: MonkeyPatchDeck) => deck._onPointerUp,
} as const;


/**
 * This is undesirable and liable to break in the future, but we need to monkeypatch
 * deck.gl to fix issues with event handling in popout windows.
 * 
 * Apply this function to a deck instance whenever the outer container changes.
 */
export function rebindMouseEvents(deckO: Deck<any>, selectionLayer?: EditableGeoJsonLayer) {
    const deck = deckO as MonkeyPatchDeck;
    //! suspected source of future problems... in order for mjolnir.js to re-bind events
    // could consider more feature-detection to alert us when deck version changes etc...
    // I think that just adds noise, we should be aware that we want to remove this & any changes
    // to deck.gl version are not to be taken lightly.

    //we tracked down the place where the event manager is created...
    //& have something based on that to re-create it.

    const oldEventManager = deck.eventManager as EventManager;
    if (!oldEventManager) {
        return;
    }
    // https://visgl.github.io/mjolnir.js/docs/api-reference/event-manager
    // > Element must be supplied when constructing EventManager and cannot be reassigned.
    // > To change the event target, destroy the existing event manager instance and construct a new one.
    // The element will be the same, but we need a new event manager so that events on window work in popout.
    // using deck.props.parent || deck.canvas; - will be the same element, but keeping copied code similar
    // const element = oldEventManager.getElement(); 
    oldEventManager.destroy();
    
    // this causes a lot of problems... but maybe if we can somewhat replicate the relevant parts of deck.ts
    // if (deck.device) deck.animationLoop.props.onInitialize(deck);

    // first make a new event manager, then make sure the viewManager knows about it
    
    // it seems like the public API of EventManager won't necessarily easily let us effectively clone the old one
    const eventManager = new EventManager(deck.props.parent || deck.canvas, {
        touchAction: deck.props.touchAction,
        recognizers: Object.keys(RECOGNIZERS).map((eventName: string) => {
            // Resolve recognizer settings
            const [RecognizerConstructor, defaultOptions, recognizeWith, requestFailure] =
                (RECOGNIZERS as any)[eventName];
            const optionsOverride = (deck.props as any).eventRecognizerOptions?.[eventName];
            const options = { ...defaultOptions, ...optionsOverride, event: eventName };
            return {
                recognizer: new RecognizerConstructor(options),
                recognizeWith,
                requestFailure
            };
        }),
        events: Object.fromEntries(
            Object.entries(DECK_EVENT_MANAGER_EVENTS)
                .map(([eventName, getHandler]) => [eventName, getHandler(deck)])
                .filter(([, handler]) => typeof handler === "function")
        ),
    });
    deck.eventManager = eventManager;
    for (const eventType in EVENT_HANDLERS) {
        deck.eventManager.on(eventType, deck._onEvent);
    }

    // we might get away with this? it uses `_eventManager` to make a `new Controller()`...
    const viewManager = deck.viewManager as any;
    viewManager._eventManager = eventManager;
    // make sure it doesn't try to re-use controllers...
    viewManager.controllers = {};
    viewManager._rebuildViewports();
    viewManager.setNeedsUpdate("MDV monkeypatch change event manager");
    if (!selectionLayer) return;
    queueMicrotask(() => {
        try {
            const activeSelectionLayer = resolveDeckSelectionLayer(deckO, selectionLayer);
            if (!activeSelectionLayer?.context?.deck) {
                return;
            }
            const layerInternal = activeSelectionLayer as {
                internalState?: unknown;
            };
            if (layerInternal.internalState) {
                return;
            }
            activeSelectionLayer.initializeState();
        } catch (e) {
            console.error("Error re-initializing (monkey-patching) editable layer state", e);
        }
    });
}