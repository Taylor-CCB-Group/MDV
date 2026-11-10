import type { DeckGLRef } from "@deck.gl/react";
import { type RefObject, useEffect } from "react";
import { rebindMouseEvents } from "@/lib/deckMonkeypatch";
import { useOuterContainer } from "../screen_state";

/**
 * Re-create deck's event manager whenever the chart moves to another window (popout / pop-in).
 *
 * mjolnir binds pointermove / pointerup to the window that owned the canvas when the event
 * manager was created, so without this a popped-out chart never sees the end of a tap or drag:
 * clicks, double-clicks and drags (including selection tools) silently do nothing.
 */
export function useRebindDeckEvents(deckRef: RefObject<DeckGLRef | null>) {
    const outerContainer = useOuterContainer();
    // biome-ignore lint/correctness/useExhaustiveDependencies: deckRef is stable; outerContainer is the trigger
    useEffect(() => {
        outerContainer; // make sure the hook runs when this changes
        const deck = deckRef.current?.deck;
        if (!deck) return;
        try {
            rebindMouseEvents(deck);
        } catch (e) {
            console.error(
                "attempt to reset deck eventManager element failed - could be related to brittle deck monkeypatch",
                e,
            );
        }
    }, [outerContainer]);
}
