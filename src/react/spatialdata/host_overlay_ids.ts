export type DeckOverlayId =
    | "grey_scatter"
    | "scatter"
    | "gate_display"
    | "selection"
    | "gate_labels";

export const DECK_OVERLAY_IDS: DeckOverlayId[] = [
    "grey_scatter",
    "scatter",
    "gate_display",
    "selection",
    "gate_labels",
];

export const DECK_OVERLAY_LABELS: Record<DeckOverlayId, string> = {
    grey_scatter: "Background scatter",
    scatter: "Scatter",
    gate_display: "Gates",
    selection: "Selection",
    gate_labels: "Gate labels",
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
