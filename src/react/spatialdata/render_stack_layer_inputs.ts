import type { RenderStack, RenderStackEntry } from "@spatialdata/layers";
import type { LayerConfig } from "@spatialdata/vis";
import { renderStackToLayerInputs, type RenderStackLayerInputs } from "@spatialdata/vis";

import { spatialEntryAsLayerConfig } from "./render_stack_entry_state";

export type RenderStackLayerInputsCache = {
    layers: Record<string, LayerConfig>;
    layerOrder: string[];
};

export function createRenderStackLayerInputsCache(): RenderStackLayerInputsCache {
    return { layers: {}, layerOrder: [] };
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
        const nextConfig = spatialEntryAsLayerConfig(
            entry as Extract<RenderStackEntry, { kind: "spatial" }>,
        );
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
