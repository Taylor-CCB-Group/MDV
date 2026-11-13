export type DeckOverlayId = "grey_scatter" | "scatter" | "gate_display" | "selection" | "gate_labels";

export const DECK_OVERLAY_IDS: DeckOverlayId[] = [
    "grey_scatter",
    "scatter",
    "gate_display",
    "selection",
    "gate_labels",
];

export const DECK_OVERLAY_LABELS: Record<DeckOverlayId, string> = {
    grey_scatter: "Filtered-out data points",
    scatter: "Data points",
    gate_display: "Gates",
    selection: "Selection",
    gate_labels: "Gate labels",
};

/**
 * Overlays the layer dialog lists. The other two follow along: the grey scatter is the
 * filtered-out half of the data points (shown or hidden with them, and switched by their
 * "Action on Filter"), and the selection overlay is driven by the selection tools.
 */
export const LISTED_DECK_OVERLAY_IDS: DeckOverlayId[] = ["scatter", "gate_display", "gate_labels"];

/** Overlays whose visibility follows another's in the layer dialog. */
export const DECK_OVERLAY_VISIBILITY_FOLLOWERS: Partial<Record<DeckOverlayId, DeckOverlayId[]>> = {
    scatter: ["grey_scatter"],
};

export function deckHostLayerId(deckId: DeckOverlayId): string {
    return `deck:${deckId}`;
}

export function isDeckHostLayerId(id: string): boolean {
    return id.startsWith("deck:");
}

export function deckIdFromHostLayerId(hostLayerId: string): DeckOverlayId | null {
    if (!hostLayerId.startsWith("deck:")) return null;
    const deckId = hostLayerId.slice("deck:".length) as DeckOverlayId;
    return DECK_OVERLAY_IDS.includes(deckId) ? deckId : null;
}
