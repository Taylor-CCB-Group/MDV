import type { RenderStack, RenderStackEntry } from "@spatialdata/layers";

import { DECK_OVERLAY_LABELS, deckIdFromHostLayerId } from "./host_overlay_ids";

export function renderStackEntryDisplayName(entry: RenderStackEntry): string {
    if (entry.kind === "host") {
        const deckId = deckIdFromHostLayerId(entry.source.hostLayerId);
        if (deckId) {
            return DECK_OVERLAY_LABELS[deckId];
        }
        return entry.source.hostLayerId;
    }
    if (entry.kind === "spatial") {
        return `${entry.source.elementType}: ${entry.source.elementKey}`;
    }
    return entry.id;
}

export function renderStackOrderLabel(stack: RenderStack): string {
    return stack.entries.map((entry) => renderStackEntryDisplayName(entry)).join(" → ");
}
