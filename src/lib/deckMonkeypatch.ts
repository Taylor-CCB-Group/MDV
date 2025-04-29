import { EventManager, InputDirection, Pan, Pinch, Tap } from 'mjolnir.js';
// import { EVENT_HANDLERS, RECOGNIZERS } from '@deck.gl/core/dist/lib/constants';
import type { Deck } from '@deck.gl/core';

type MonkeyPatchDeck = Deck<any> & {
    eventManager: EventManager;
    viewManager: any;
    canvas: HTMLCanvasElement;
}

const EVENT_HANDLERS: { [eventName: string]: string } = {
    click: 'onClick',
    panstart: 'onDragStart',
    panmove: 'onDrag',
    panend: 'onDragEnd'
} as const;

const RECOGNIZERS = {
    multipan: [Pan, { threshold: 10, direction: InputDirection.Vertical, pointers: 2 }],
    pinch: [Pinch, {}, null, ['multipan']],
    pan: [Pan, { threshold: 1 }, ['pinch'], ['multipan']],
    dblclick: [Tap, { event: 'dblclick', taps: 2 }],
    click: [Tap, { event: 'click' }, null, ['dblclick']]
} as const;



/**
 * This is undesirable and liable to break in the future, but we need to monkeypatch
 * deck.gl to fix issues with event handling in popout windows.
 * 
 * Apply this function to a deck instance whenever the outer container changes.
 */
export function rebindMouseEvents(deckO: Deck<any>) {
    const deck = deckO as MonkeyPatchDeck;
    //! suspected source of future problems... in order for mjolnir.js to re-bind events
    //we tracked down the place where the event manager is created...
    //but this is frought with problems
    //- EditableLayer modes not working
    //- lots of failed assertions in deck.gl
    //- glitchiness from a user perspective with pan/zoom etc after switching/back...

    // we pass deck not because it's the right type, but because we know it will just look at the device property.
    // deck.viewManager?.setNeedsUpdate("MDV useOuterContainer() changed (fullscreen/popout)");
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
        events: {
            pointerdown: deck._onPointerDown,
            pointermove: deck._onPointerMove,
            pointerleave: deck._onPointerMove
        }
    });
    for (const eventType in EVENT_HANDLERS) {
        deck.eventManager.on(eventType, deck._onEvent);
    }

    // we might get away with this? it uses `_eventManager` to make a `new Controller()` in 
    deck.viewManager._eventManager = eventManager;
    deck.viewManager.controllers = {};
    deck.viewManager._rebuildViewports();
    deck.viewManager.setNeedsUpdate("MDV monkeypatch change event manager");
}