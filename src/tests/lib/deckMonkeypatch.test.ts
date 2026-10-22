import { Deck } from "@deck.gl/core";
import { EditableGeoJsonLayer } from "@deck.gl-community/editable-layers";
import { describe, expect, it, vi } from "vitest";
import { rebindMouseEvents } from "@/lib/deckMonkeypatch";

const DEFAULT_CANVAS_ID = "default-canvas";

/** A Deck with its event-related state stubbed, as it is once its device has initialized. */
function makeDeck() {
    const oldEventManager = { destroy: vi.fn(), on: vi.fn() };
    const newEventManager = { destroy: vi.fn(), on: vi.fn() };
    const canvas = document.createElement("canvas");
    const eventRoot = document.createElement("div");
    const viewManager = { _eventManager: oldEventManager, setProps: vi.fn() };
    const layers: unknown[] = [];
    const deck = Object.assign(Object.create(Deck.prototype), {
        canvas,
        eventManager: oldEventManager,
        eventManagers: { [DEFAULT_CANVAS_ID]: oldEventManager },
        viewManager,
        layerManager: { getLayers: () => layers },
        _isMultiCanvasMode: () => false,
        _getEventRoot: vi.fn(() => eventRoot),
        _createEventManager: vi.fn(() => newEventManager),
    });
    return { deck, oldEventManager, newEventManager, canvas, eventRoot, viewManager, layers };
}

describe("deckMonkeypatch", () => {
    it("relies on deck internals that still exist", () => {
        // if a deck update renames these, rebindMouseEvents bails and popped-out charts lose mouse events
        for (const name of ["_createEventManager", "_getEventRoot", "_isMultiCanvasMode"]) {
            expect(typeof Reflect.get(Deck.prototype, name)).toBe("function");
        }
    });

    it("creates the replacement manager with deck's own configuration and event root", () => {
        const { deck, canvas, eventRoot } = makeDeck();

        rebindMouseEvents(deck);

        expect(deck._getEventRoot).toHaveBeenCalledWith(canvas);
        expect(deck._createEventManager).toHaveBeenCalledWith(eventRoot);
    });

    it("replaces the manager everywhere deck holds it, so controllers are rebuilt on it", () => {
        const { deck, oldEventManager, newEventManager, viewManager } = makeDeck();

        rebindMouseEvents(deck);

        expect(deck.eventManager).toBe(newEventManager);
        expect(deck.eventManagers[DEFAULT_CANVAS_ID]).toBe(newEventManager);
        expect(viewManager._eventManager).toBe(newEventManager);
        expect(viewManager.setProps).toHaveBeenCalledWith({ eventManagers: { [DEFAULT_CANVAS_ID]: newEventManager } });
        // old controllers unsubscribe from the old manager while the view manager rebuilds
        expect(oldEventManager.destroy).toHaveBeenCalledOnce();
        expect(oldEventManager.destroy.mock.invocationCallOrder[0]).toBeGreaterThan(
            viewManager.setProps.mock.invocationCallOrder[0],
        );
    });

    it("re-attaches editable layers found in the deck to the new manager", () => {
        // the instance in the deck may be a clone of the one the chart created (e.g. SpatialData.js
        // re-ids host layers), so it has to be found in the deck rather than passed in
        const { deck, newEventManager, layers } = makeDeck();
        const eventHandler = vi.fn();
        const editable = new EditableGeoJsonLayer({ id: "deck:selection-#detail#", data: { type: "FeatureCollection", features: [] } });
        Object.assign(editable, { context: { deck }, state: { _editableLayerState: { eventHandler } } });
        layers.push(editable);

        rebindMouseEvents(deck);

        const types = newEventManager.on.mock.calls.filter(([, handler]) => handler === eventHandler).map(([type]) => type);
        expect(types).toEqual(expect.arrayContaining(["click", "dblclick", "pointermove", "panstart", "panmove", "panend"]));
        expect(new Set(types).size).toBe(types.length);
    });

    it("does nothing before deck has created its event manager", () => {
        const { deck } = makeDeck();
        deck.eventManager = null;

        rebindMouseEvents(deck);

        expect(deck._createEventManager).not.toHaveBeenCalled();
    });
});
