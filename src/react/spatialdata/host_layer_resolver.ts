import type { RenderStack, RenderStackHostEntry } from "@spatialdata/layers";
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

export function renderStackHostFingerprint(stack: RenderStack | undefined): string {
    if (!stack) return "";
    return stack.entries
        .filter((entry): entry is RenderStackHostEntry => entry.kind === "host")
        .map((entry) => `${entry.id}:${entry.visible}`)
        .join("|");
}

const hostLayerCloneCache = new WeakMap<Layer, Map<string, Layer>>();

function cloneHostLayer(source: Layer, entryId: string): Layer {
    let clonesForSource = hostLayerCloneCache.get(source);
    if (!clonesForSource) {
        clonesForSource = new Map();
        hostLayerCloneCache.set(source, clonesForSource);
    }
    let clone = clonesForSource.get(entryId);
    if (!clone) {
        clone = source.clone({ id: entryId }) as Layer;
        clonesForSource.set(entryId, clone);
    }
    return clone;
}

/**
 * Resolve visible host stack entries to deck layers, cloning each source layer once per
 * entry id. Cosmetic spatial-layer edits (e.g. opacity) can then refresh the render
 * stack shell without re-cloning scatter/gate overlays on every frame.
 */
export function resolveCachedHostDeckLayers(
    stack: RenderStack | undefined,
    resolver: (entry: RenderStackHostEntry) => Layer | Layer[] | null | undefined,
): Layer[] {
    if (!stack) return [];

    const layers: Layer[] = [];
    for (const entry of stack.entries) {
        if (entry.kind !== "host" || !entry.visible) continue;
        const resolved = resolver(entry);
        if (!resolved) continue;
        const sources = Array.isArray(resolved) ? resolved : [resolved];
        for (const source of sources) {
            if (!source) continue;
            layers.push(cloneHostLayer(source, entry.id));
        }
    }
    return layers;
}
