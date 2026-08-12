import type { Layer } from "@deck.gl/core";
import type {
    RenderStack,
    RenderStackEntry,
    RenderStackHostEntry,
} from "@spatialdata/layers";
import {
    layerConfig,
    renderStackOrder,
    type LayerConfig,
    type RenderStackLayerInputs,
} from "@spatialdata/vis";
import { useMemo, useRef } from "react";

import { withoutMdvFieldSpecs } from "./field_spec_projection";
import { deckIdFromHostLayerId, type DeckOverlayId } from "./host_overlay_ids";
import { touchRenderStack, renderStackSpatialRevision } from "./render_stack_observe";
import { measureSpatial } from "./perf";

export type MdvDeckOverlayLayers = Record<DeckOverlayId, Layer | null>;

export type RenderStackLayerInputsCache = {
    layers: Record<string, LayerConfig>;
    layerOrder: string[];
    layerConfigSignatures: Record<string, string>;
};

export function createRenderStackLayerInputsCache(): RenderStackLayerInputsCache {
    return { layers: {}, layerOrder: [], layerConfigSignatures: {} };
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

function stableSignatureValue(value: unknown): unknown {
    if (Array.isArray(value)) return value.map(stableSignatureValue);
    if (!value || typeof value !== "object") return value;

    const record = value as Record<string, unknown>;
    return Object.fromEntries(
        Object.keys(record)
            .sort()
            .map((key) => [key, stableSignatureValue(record[key])]),
    );
}

function layerConfigReplacementSignature(config: LayerConfig): string {
    const fillColorByColumn =
        "fillColorByColumn" in config ? config.fillColorByColumn : undefined;
    const tooltipFields = "tooltipFields" in config ? config.tooltipFields : undefined;
    return JSON.stringify({
        type: config.type,
        elementKey: config.elementKey,
        fillColorByColumn: stableSignatureValue(fillColorByColumn),
        tooltipFields: stableSignatureValue(tooltipFields),
    });
}

/**
 * Keep a stable `layers` object identity across cosmetic edits so
 * `useLayerData` does not re-enter async geometry loads. Mutate layer configs
 * in place for cosmetic props; replace the individual layer config when
 * table/annotation-driving props change so downstream projection caches see a
 * clear A -> B transition.
 */
export function syncRenderStackLayerInputs(
    stack: RenderStack,
    cache: RenderStackLayerInputsCache,
): RenderStackLayerInputs {
    const nextIds = new Set<string>();

    for (const entry of stack.entries) {
        if (entry.kind !== "spatial") continue;
        nextIds.add(entry.id);
        // MDV's field specs stop here: past this point a layer config is the
        // viewer's, and a `RowsAsColsQuery` is not something to hand it.
        const nextConfig = withoutMdvFieldSpecs(spatialEntryAsLayerConfig(entry));
        const nextSignature = layerConfigReplacementSignature(nextConfig);
        const existing = cache.layers[entry.id];
        if (existing && cache.layerConfigSignatures[entry.id] === nextSignature) {
            Object.assign(existing, nextConfig);
        } else {
            cache.layers[entry.id] = nextConfig;
            cache.layerConfigSignatures[entry.id] = nextSignature;
        }
    }

    for (const id of Object.keys(cache.layers)) {
        if (!nextIds.has(id)) {
            delete cache.layers[id];
            delete cache.layerConfigSignatures[id];
        }
    }

    const nextOrder = renderStackOrder(stack, cache.layerOrder);
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

function cloneHostLayer(source: Layer, entryId: string): Layer {
    return source.clone({ id: entryId });
}

/**
 * Resolve visible host stack entries to deck layers. The source layers are owned by
 * MDV hooks, so clone at each adapter refresh to preserve the latest source props.
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
    propsGeneration,
    hostLayerResolver,
}: {
    stack: RenderStack | undefined;
    generation: number;
    propsGeneration: number;
    hostLayerResolver: ReturnType<typeof createMdvHostLayerResolver>;
}) {
    const layerInputsCacheRef = useRef(createRenderStackLayerInputsCache());
    // `measureSpatial` is a no-op unless localStorage.MDV_SPATIAL_PERF === "1".
    // The `count` of these labels == number of adapter renders during a capture,
    // i.e. how often a cosmetic image edit re-renders SpatialDataViewer.
    measureSpatial("adapter.observe", () => observeRenderStack(stack));
    const spatialRevision = measureSpatial("adapter.revision", () =>
        `${propsGeneration}:${renderStackSpatialRevision(stack)}`,
    );

    const synced = measureSpatial("adapter.sync", () =>
        !stack
            ? { layers: {}, layerOrder: [] as string[] }
            : syncRenderStackLayerInputs(stack, layerInputsCacheRef.current),
    );

    // useLayerData (upstream) memoizes Viv props on the layers *record* reference.
    // In-place config mutation keeps that reference stable, so shallow-copy the map
    // when spatial revision changes to invalidate Viv without replacing layer configs.
    const layers = useMemo(() => {
        void spatialRevision; // trigger-only dep: re-shallow-copy to invalidate the upstream Viv memo
        return { ...synced.layers };
    }, [spatialRevision, synced.layers]);

    const hostFingerprint = renderStackHostFingerprint(stack);
    const deckLayers = useMemo(() => {
        void hostFingerprint; // trigger-only dep: re-clone host layers when host visibility changes
        void generation; // trigger-only dep: re-clone on coarse render-stack generation bumps
        return resolveCachedHostDeckLayers(stack, hostLayerResolver);
    }, [stack, hostFingerprint, hostLayerResolver, generation]);

    return {
        layers,
        layerOrder: synced.layerOrder,
        deckLayers,
    };
}
