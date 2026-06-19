import type { RenderStackHostEntry } from "@spatialdata/layers";
import type { Layer } from "@deck.gl/core";

import { deckIdFromHostLayerId, type DeckOverlayId } from "./host_overlay_ids";

export type MdvDeckOverlayLayers = Record<DeckOverlayId, Layer | null>;

export function createMdvHostLayerResolver(overlays: MdvDeckOverlayLayers) {
    return (entry: RenderStackHostEntry): Layer | null => {
        const deckId = deckIdFromHostLayerId(entry.source.hostLayerId);
        if (!deckId) return null;
        return overlays[deckId] ?? null;
    };
}
