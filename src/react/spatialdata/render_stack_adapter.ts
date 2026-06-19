import type { Layer } from "@deck.gl/core";
import type {
    RenderStack,
    RenderStackEntry,
    RenderStackHostEntry,
} from "@spatialdata/layers";
import {
    layerConfig,
    renderStackToLayerInputs,
    type LayerConfig,
    type RenderStackLayerInputs,
} from "@spatialdata/vis";
import { useMemo, useRef } from "react";

import { deckIdFromHostLayerId, type DeckOverlayId } from "./host_overlay_ids";

export type MdvDeckOverlayLayers = Record<DeckOverlayId, Layer | null>;

export type RenderStackLayerInputsCache = {
    layers: Record<string, LayerConfig>;
    layerOrder: string[];
};

export function createRenderStackLayerInputsCache(): RenderStackLayerInputsCache {
    return { layers: {}, layerOrder: [] };
}

export function createMdvHostLayerResolver(overlays: MdvDeckOverlayLayers) {
    return (entry: RenderStackHostEntry): Layer | null => {
        const deckId = deckIdFromHostLayerId(entry.source.hostLayerId);
        if (!deckId) return null;
        return overlays[deckId] ?? null;
    };
}

export function spatialEntryAsLayerConfig(
    entry: Extract<RenderStackEntry, { kind: "spatial" }>,
): LayerConfig {
    return layerConfig(
        entry.source.elementType,
        {
            id: entry.id,
            elementKey: entry.source.elementKey,
            visible: entry.visible,
            opacity: typeof entry.props.opacity === "number" ? entry.props.opacity : 1,
        },
        entry.props,
    );
}

/**
 * Keep a stable `layers` object identity across cosmetic edits so
 * `useLayerData` does not re-enter async geometry loads. Mutate layer configs
 * in place when props change; only replace/add/remove entries structurally.
 */
export function syncRenderStackLayerInputs(
    stack: RenderStack,
    cache: RenderStackLayerInputsCache,
): RenderStackLayerInputs {
    const nextIds = new Set<string>();

    for (const entry of stack.entries) {
        if (entry.kind !== "spatial") continue;
        nextIds.add(entry.id);
        const nextConfig = spatialEntryAsLayerConfig(entry);
        const existing = cache.layers[entry.id];
        if (existing) {
            Object.assign(existing, nextConfig);
        } else {
            cache.layers[entry.id] = nextConfig;
        }
    }

    for (const id of Object.keys(cache.layers)) {
        if (!nextIds.has(id)) {
            delete cache.layers[id];
        }
    }

    const nextOrder = renderStackToLayerInputs(stack).layerOrder;
    if (nextOrder.join("\0") !== cache.layerOrder.join("\0")) {
        cache.layerOrder = nextOrder;
    }

    return { layers: cache.layers, layerOrder: cache.layerOrder };
}

function renderStackHostFingerprint(stack: RenderStack | undefined): string {
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
 * entry id. Cosmetic spatial-layer edits (e.g. opacity) can then refresh adapter outputs
 * without re-cloning scatter/gate overlays on every frame.
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

function observeSpatialRenderStack(stack: RenderStack | undefined) {
    for (const entry of stack?.entries ?? []) {
        if (entry.kind !== "spatial") continue;
        void entry.id;
        void entry.visible;
        void entry.props;
    }
}

export function useRenderStackAdapter({
    stack,
    generation,
    hostLayerResolver,
}: {
    stack: RenderStack | undefined;
    generation: number;
    hostLayerResolver: ReturnType<typeof createMdvHostLayerResolver>;
}) {
    const layerInputsCacheRef = useRef(createRenderStackLayerInputsCache());
    observeSpatialRenderStack(stack);

    const layerInputs = useMemo(() => {
        if (!stack) return { layers: {}, layerOrder: [] as string[] };
        return syncRenderStackLayerInputs(stack, layerInputsCacheRef.current);
    }, [stack, generation]);

    const hostFingerprint = renderStackHostFingerprint(stack);
    const deckLayers = useMemo(
        () => resolveCachedHostDeckLayers(stack, hostLayerResolver),
        [stack, hostFingerprint, hostLayerResolver],
    );

    return {
        layers: layerInputs.layers,
        layerOrder: layerInputs.layerOrder,
        deckLayers,
    };
}
