import { useRef } from "react";

/**
 * Increment when deck layer ownership mode toggles (overlay ↔ grid).
 * deck.gl finalizes removed layers (`assert(!this.internalState)`); clones must be
 * fresh JS instances after each toggle, not reused from a prior mount.
 */
export function useDeckMountEpoch(modeActive: boolean): number {
    const epochRef = useRef(0);
    const prevRef = useRef(modeActive);
    if (prevRef.current !== modeActive) {
        epochRef.current += 1;
        prevRef.current = modeActive;
    }
    return epochRef.current;
}
