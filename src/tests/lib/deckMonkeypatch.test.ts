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

import { rebindMouseEvents, resolveDeckSelectionLayer } from "@/lib/deckMonkeypatch";

describe("deckMonkeypatch", () => {
    it("registers pointerup when rebinding deck events", () => {
        const deck = {
            props: { parent: null, touchAction: "none", eventRecognizerOptions: {} },
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

        rebindMouseEvents(deck);

        const ctorArgs = eventManagerCtor.mock.calls.at(-1);
        expect(ctorArgs).toBeTruthy();
        const options = ctorArgs?.[1] as { events: Record<string, unknown> };
        expect(options.events.pointerdown).toBe(deck._onPointerDown);
        expect(options.events.pointermove).toBe(deck._onPointerMove);
        expect(options.events.pointerleave).toBe(deck._onPointerMove);
        expect(options.events.pointerup).toBe(deck._onPointerUp);
    });

    it("resolves mounted selection layers before re-initializing editable state", () => {
        const mountedSelection = {
            id: "selection_-#chartchart-array-grid#",
            context: { deck: {} },
            initializeState: vi.fn(),
        };
        const staleSelection = {
            id: "selection_-#chartdetail-react#",
            initializeState: vi.fn(),
        };
        const deck = {
            layerManager: {
                layers: [mountedSelection],
            },
        } as any;

        expect(resolveDeckSelectionLayer(deck, staleSelection as any)).toBe(mountedSelection);
    });
});
