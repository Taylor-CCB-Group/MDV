import type { Layer } from "@deck.gl/core";
import {
    type RenderStack,
    type RenderStackEntry,
    type RenderStackHostEntry,
    resolveRenderStackGroupProps,
    walkRenderStackEntries,
} from "@spatialdata/layers";
import {
    type LayerConfig,
    layerConfig,
    type RenderStackLayerInputs,
    renderStackOrder,
    renderStackToLayerInputs,
    resolveRenderStackHostLayers,
} from "@spatialdata/vis";
import { useMemo, useRef } from "react";

import { withoutMdvFieldSpecs } from "./field_spec_projection";
import { type DeckOverlayId, deckIdFromHostLayerId } from "./host_overlay_ids";
import { measureSpatial } from "./perf";
import { DISPLAY_NAME_PROP } from "./render_stack_display";
import { renderStackSpatialRevision, touchRenderStack } from "./render_stack_observe";

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

export function spatialEntryAsLayerConfig(entry: Extract<RenderStackEntry, { kind: "spatial" }>): LayerConfig {
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
 * Take MDV's own bookkeeping off a layer config on its way to the viewer: the field
 * specs (see `withoutMdvFieldSpecs`) and the user's name for the layer, which only the
 * layer dialog reads.
 */
export function withoutMdvOnlyProps<T extends object>(config: T): T {
    const stripped = withoutMdvFieldSpecs(config);
    if (!(DISPLAY_NAME_PROP in stripped)) return stripped;
    const { [DISPLAY_NAME_PROP]: _removed, ...rest } = stripped as T & { [DISPLAY_NAME_PROP]?: unknown };
    return rest as T;
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
    const fillColorByColumn = "fillColorByColumn" in config ? config.fillColorByColumn : undefined;
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
    // Spatial entries at any depth, with their groups' visibility and opacity applied.
    const inputs = renderStackToLayerInputs(stack);

    for (const id of inputs.layerOrder) {
        nextIds.add(id);
        // MDV's field specs and layer names stop here: past this point a layer config
        // is the viewer's, and a `RowsAsColsQuery` is not something to hand it.
        const nextConfig = withoutMdvOnlyProps(inputs.layers[id]);
        const nextSignature = layerConfigReplacementSignature(nextConfig);
        const existing = cache.layers[id];
        if (existing && cache.layerConfigSignatures[id] === nextSignature) {
            Object.assign(existing, nextConfig);
        } else {
            cache.layers[id] = nextConfig;
            cache.layerConfigSignatures[id] = nextSignature;
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

/** What host layers depend on besides their sources: entry visibility, and the
 *  visibility and opacity of the groups above them. */
function renderStackHostFingerprint(stack: RenderStack | undefined): string {
    if (!stack) return "";
    const parts: string[] = [];
    for (const entry of walkRenderStackEntries(stack.entries)) {
        if (entry.kind === "host") parts.push(`${entry.id}:${entry.visible}`);
        if (entry.kind === "group") {
            const { opacity, isolate } = resolveRenderStackGroupProps(entry);
            parts.push(`${entry.id}:${entry.visible}:${opacity}:${isolate}:${entry.children.length}`);
        }
    }
    return parts.join("|");
}

/**
 * Resolve visible host stack entries, at any depth, to deck layers. The source layers
 * are owned by MDV hooks, so clone at each adapter refresh to preserve the latest
 * source props.
 */
export function resolveCachedHostDeckLayers(
    stack: RenderStack | undefined,
    resolver: (entry: RenderStackHostEntry) => Layer | Layer[] | null | undefined,
): Layer[] {
    return resolveRenderStackHostLayers(stack, resolver);
}

function entryStructure(entry: RenderStackEntry): RenderStackEntry {
    if (entry.kind !== "group") return { ...entry, props: {} };
    const { opacity, blendMode, isolate } = entry.props;
    return {
        kind: "group",
        id: entry.id,
        visible: entry.visible,
        props: { opacity, blendMode, isolate },
        children: entry.children.map(entryStructure),
    };
}

/**
 * The stack's tree for the viewer to place layers by: ids, sources, visibility and
 * group props, with entry props left out. The observable stack is edited in place, so
 * this is a plain copy, and keeps its identity until its value changes: a cosmetic
 * edit to an entry's props must not hand the viewer a new stack. Call from an
 * observer, so structural edits re-render.
 */
function renderStackStructure(stack: RenderStack): RenderStack {
    return { schemaVersion: stack.schemaVersion, entries: stack.entries.map(entryStructure) };
}

export function useRenderStackStructure(stack: RenderStack | undefined): RenderStack | undefined {
    // Building the key reads the whole structure, so an observer re-renders on any
    // structural edit; the copy is rebuilt only when that key changes.
    const key = stack ? JSON.stringify(renderStackStructure(stack)) : "";
    return useMemo(() => {
        void key; // trigger-only dep: the stack is mutated in place, so its identity never changes
        return stack && renderStackStructure(stack);
    }, [stack, key]);
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
    // `measureSpatial` is a no-op unless localStorage.MDV_SPATIAL_PERF === "1" —
    // see `perf.ts` for the labels, what a healthy capture looks like, and the
    // warning about leaving that flag set.
    measureSpatial("adapter.observe", () => observeRenderStack(stack));
    const spatialRevision = measureSpatial(
        "adapter.revision",
        () => `${propsGeneration}:${renderStackSpatialRevision(stack)}`,
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
