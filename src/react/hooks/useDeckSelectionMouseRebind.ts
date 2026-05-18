import { useEffect, type RefObject } from "react";
import type { Deck } from "@deck.gl/core";
import type { EditableGeoJsonLayer } from "@deck.gl-community/editable-layers";
import { rebindMouseEvents } from "@/lib/deckMonkeypatch";

type DeckSelectionMouseRebindOptions = {
    /** When false, skip binding (e.g. overlay Deck unmounted while grid is shown). */
    enabled?: boolean;
    /** Bust the effect when the Deck canvas remounts (overlay vs density grid). */
    canvasKey?: string;
};

type DeckRefTarget = {
    deck?: Deck;
} | null;

const REBIND_POLL_MS = 50;
const REBIND_POLL_MAX_MS = 3000;

/**
 * Rebind deck + editable-layer mouse handlers when the chart moves to a popout
 * or when switching between single-view and density-grid Deck canvases.
 */
export function useDeckSelectionMouseRebind(
    outerContainer: HTMLElement | undefined,
    selectionLayer: EditableGeoJsonLayer | undefined,
    deckRef: RefObject<DeckRefTarget>,
    options?: DeckSelectionMouseRebindOptions,
) {
    const enabled = options?.enabled !== false;
    const canvasKey = options?.canvasKey ?? "default";

    useEffect(() => {
        if (!enabled) return;
        outerContainer;

        let cancelled = false;
        let cleanup: (() => void) | undefined;

        const tryBind = () => {
            const deck = deckRef.current?.deck;
            if (!deck || cancelled) return false;
            try {
                cleanup = rebindMouseEvents(deck, selectionLayer);
                return true;
            } catch (e) {
                console.error(
                    "attempt to reset deck eventManager element failed - could be related to brittle deck monkeypatch",
                    e,
                );
                return false;
            }
        };

        if (tryBind()) {
            return () => {
                cancelled = true;
                cleanup?.();
            };
        }

        const started = Date.now();
        const timer = window.setInterval(() => {
            if (tryBind() || Date.now() - started >= REBIND_POLL_MAX_MS) {
                window.clearInterval(timer);
            }
        }, REBIND_POLL_MS);

        return () => {
            cancelled = true;
            window.clearInterval(timer);
            cleanup?.();
        };
    }, [outerContainer, selectionLayer, deckRef, enabled, canvasKey]);
}
