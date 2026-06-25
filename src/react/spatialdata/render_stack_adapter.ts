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
import { touchRenderStack, renderStackSpatialRevision } from "./render_stack_observe";
import { measureSpatial } from "./perf";

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

function observeRenderStack(stack: RenderStack | undefined) {
    touchRenderStack(stack);
    void renderStackSpatialRevision(stack);
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
    // `measureSpatial` is a no-op unless localStorage.MDV_SPATIAL_PERF === "1".
    // The `count` of these labels == number of adapter renders during a capture,
    // i.e. how often a cosmetic image edit re-renders SpatialDataViewer.
    measureSpatial("adapter.observe", () => observeRenderStack(stack));
    const spatialRevision = measureSpatial("adapter.revision", () =>
        renderStackSpatialRevision(stack),
    );
    void generation;

    const synced = measureSpatial("adapter.sync", () =>
        !stack
            ? { layers: {}, layerOrder: [] as string[] }
            : syncRenderStackLayerInputs(stack, layerInputsCacheRef.current),
    );

    // useLayerData (upstream) memoizes Viv props on the layers *record* reference.
    // In-place config mutation keeps that reference stable, so shallow-copy the map
    // when spatial revision changes to invalidate Viv without replacing layer configs.
    const layers = useMemo(
        () => ({ ...synced.layers }),
        [spatialRevision, synced.layers],
    );

    const hostFingerprint = renderStackHostFingerprint(stack);
    const deckLayers = useMemo(
        () => resolveCachedHostDeckLayers(stack, hostLayerResolver),
        [stack, hostFingerprint, hostLayerResolver, generation],
    );

    return {
        layers,
        layerOrder: synced.layerOrder,
        deckLayers,
    };
}
