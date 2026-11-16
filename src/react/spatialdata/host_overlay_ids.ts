export type DeckOverlayId = "grey_scatter" | "scatter" | "density" | "gate_display" | "selection" | "gate_labels";

/** In drawing order, bottom first. A saved stack missing one gets it in this place. */
export const DECK_OVERLAY_IDS: DeckOverlayId[] = [
    "grey_scatter",
    "scatter",
    "density",
    "gate_display",
    "selection",
    "gate_labels",
];

export const DECK_OVERLAY_LABELS: Record<DeckOverlayId, string> = {
    grey_scatter: "Filtered-out data points",
    scatter: "Data points",
    density: "Density",
    gate_display: "Gate outlines",
    selection: "Selection",
    gate_labels: "Gate labels",
};

/** A row of the layer dialog's "MDV data" list: one overlay, or a named group of them. */
export type ListedDeckOverlay =
    | { kind: "overlay"; deckId: DeckOverlayId }
    | { kind: "group"; label: string; children: DeckOverlayId[] };

/**
 * What the layer dialog lists, in order. The other two overlays follow along: the grey
 * scatter is the filtered-out half of the data points (shown or hidden with them, and
 * switched by their "Action on Filter"), and the selection overlay is driven by the
 * selection tools. A group exists only in the dialog; the stack keeps its overlays flat.
 */
export const LISTED_DECK_OVERLAYS: ListedDeckOverlay[] = [
    { kind: "overlay", deckId: "scatter" },
    { kind: "overlay", deckId: "density" },
    {
        kind: "group",
        label: "Gates",
        children: ["gate_display", "gate_labels"],
    },
];

export const LISTED_DECK_OVERLAY_IDS: DeckOverlayId[] = LISTED_DECK_OVERLAYS.flatMap((row) =>
    row.kind === "overlay" ? [row.deckId] : row.children,
);

/** Overlays whose visibility follows another's in the layer dialog. */
export const DECK_OVERLAY_VISIBILITY_FOLLOWERS: Partial<Record<DeckOverlayId, DeckOverlayId[]>> = {
    scatter: ["grey_scatter"],
};

/**
 * The group the MDV data points (grey and coloured scatter) are composited in, while their
 * compositing is not the default. Not a user group: the dialog lists it as "Data points".
 */
export const DATA_POINTS_GROUP_ID = "group:data-points";

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
