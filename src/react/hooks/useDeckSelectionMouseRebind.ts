import { rebindMouseEvents } from "@/lib/deckMonkeypatch";
import type { EditableGeoJsonLayer } from "@deck.gl-community/editable-layers";
import type { Deck } from "@deck.gl/core";
import { type RefObject, useEffect, useRef } from "react";
import { useRange } from "../spatial_context";

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
 *
 * The selection layer instance is recreated when geojson data/mode props change;
 * that must not trigger a rebind (it destroys EventManager mid-gesture). Keep the
 * latest layer in a ref and only rebind for container/canvas/enabled changes.
 */
export function useDeckSelectionMouseRebind(
    outerContainer: HTMLElement | undefined,
    selectionLayer: EditableGeoJsonLayer | undefined,
    deckRef: RefObject<DeckRefTarget>,
    options?: DeckSelectionMouseRebindOptions,
) {
    const enabled = options?.enabled !== false;
    const canvasKey = options?.canvasKey ?? "default";
    const { editingGateId } = useRange();
    const selectionLayerRef = useRef(selectionLayer);
    selectionLayerRef.current = selectionLayer;

    // biome-ignore lint/correctness/useExhaustiveDependencies: selectionLayer is recreated on every geojson edit; listing it rebinds mid-gesture (use selectionLayerRef). outerContainer/canvasKey bust the effect when the Deck host changes.
    useEffect(() => {
        if (!enabled) return;
        if (editingGateId) return;

        let cancelled = false;
        let cleanup: (() => void) | undefined;

        const tryBind = () => {
            const deck = deckRef.current?.deck;
            const layer = selectionLayerRef.current;
            if (!deck || cancelled) return false;
            try {
                cleanup = rebindMouseEvents(deck, layer);
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
    }, [outerContainer, deckRef, enabled, canvasKey, editingGateId]);
}
