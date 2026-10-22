import { describe, expect, it, vi } from "vitest";

const {
    eventManagerCtor,
    eventManagerDestroy,
    eventManagerOn,
    recognizerCtor,
} = vi.hoisted(() => {
    return {
        eventManagerCtor: vi.fn(),
        eventManagerDestroy: vi.fn(),
        eventManagerOn: vi.fn(),
        recognizerCtor: vi.fn(function recognizer(this: object) {
            return this;
        }),
    };
});

vi.mock("mjolnir.js", () => {
    class MockEventManager {
        constructor(...args: unknown[]) {
            eventManagerCtor(...args);
        }
        destroy() {
            eventManagerDestroy();
        }
        on(...args: unknown[]) {
            eventManagerOn(...args);
        }
    }
    return {
        EventManager: MockEventManager,
        InputDirection: { Vertical: "vertical" },
        Pan: recognizerCtor,
        Pinch: recognizerCtor,
        Tap: recognizerCtor,
    };
});

import { EditableGeoJsonLayer } from "@deck.gl-community/editable-layers";
import { rebindMouseEvents } from "@/lib/deckMonkeypatch";

function makeDeck(parent: HTMLElement | null = null) {
    return {
        props: { parent, touchAction: "none", eventRecognizerOptions: {} },
        canvas: {} as HTMLElement,
        eventManager: {
            destroy: vi.fn(),
        },
        viewManager: {
            _eventManager: null,
            controllers: {},
            _rebuildViewports: vi.fn(),
            setNeedsUpdate: vi.fn(),
        },
        _onPointerDown: vi.fn(),
        _onPointerMove: vi.fn(),
        _onPointerUp: vi.fn(),
        _onEvent: vi.fn(),
    } as any;
}

type RecognizerItem = { recognizeWith?: string[]; requireFailure?: string[] };

describe("deckMonkeypatch", () => {
    it("registers pointerup when rebinding deck events", () => {
        const deck = makeDeck();

        rebindMouseEvents(deck);

        const ctorArgs = eventManagerCtor.mock.calls.at(-1);
        expect(ctorArgs).toBeTruthy();
        const options = ctorArgs?.[1] as { events: Record<string, unknown> };
        expect(options.events.pointerdown).toBe(deck._onPointerDown);
        expect(options.events.pointermove).toBe(deck._onPointerMove);
        expect(options.events.pointerleave).toBe(deck._onPointerMove);
        expect(options.events.pointerup).toBe(deck._onPointerUp);
    });

    it("makes click wait for dblclick to fail, as deck does", () => {
        // mjolnir reads `requireFailure`; otherwise a double-click also emits two clicks
        rebindMouseEvents(makeDeck());

        const options = eventManagerCtor.mock.calls.at(-1)?.[1] as { recognizers: RecognizerItem[] };
        const [multipan, pinch, pan, dblclick, click] = options.recognizers;
        expect(click.requireFailure).toEqual(["dblclick"]);
        expect(pan.requireFailure).toEqual(["multipan"]);
        expect(pinch.requireFailure).toEqual(["multipan"]);
        expect(multipan.requireFailure).toBeFalsy();
        expect(dblclick.requireFailure).toBeFalsy();
        expect(options.recognizers.every((r) => !("requestFailure" in r))).toBe(true);
    });

    it("forwards click and dblclick to deck's onClick handling", () => {
        const deck = makeDeck();
        rebindMouseEvents(deck);

        const registered = eventManagerOn.mock.calls.filter(([, handler]) => handler === deck._onEvent).map(([type]) => type);
        expect(registered).toEqual(expect.arrayContaining(["click", "dblclick", "panstart", "panmove", "panend"]));
    });

    it("binds to deck's event root inside the parent, as deck does", () => {
        const parent = document.createElement("div");
        const eventRoot = document.createElement("div");
        eventRoot.className = "deck-events-root";
        parent.appendChild(eventRoot);

        rebindMouseEvents(makeDeck(parent));

        expect(eventManagerCtor.mock.calls.at(-1)?.[0]).toBe(eventRoot);
    });

    it("re-attaches editable layers found in the deck to the new event manager", () => {
        // the instance in the deck may be a clone of the one the chart created (e.g. SpatialData.js
        // re-ids host layers), so it has to be found in the deck rather than passed in
        const deck = makeDeck();
        const eventHandler = vi.fn();
        const editable = new EditableGeoJsonLayer({ id: "deck:selection-#detail#", data: { type: "FeatureCollection", features: [] } });
        Object.assign(editable, { context: { deck }, state: { _editableLayerState: { eventHandler } } });
        deck.layerManager = { getLayers: () => [editable] };
        eventManagerOn.mockClear();

        rebindMouseEvents(deck);

        const types = eventManagerOn.mock.calls.filter(([, handler]) => handler === eventHandler).map(([type]) => type);
        expect(types).toEqual(expect.arrayContaining(["click", "dblclick", "pointermove", "panstart", "panmove", "panend"]));
        expect(new Set(types).size).toBe(types.length);
    });
});
