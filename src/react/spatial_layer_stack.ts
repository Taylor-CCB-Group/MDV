import type { SpatialData } from "@spatialdata/core";
import type {
    ElementsByType,
    LayerConfig,
    LayerType,
} from "@spatialdata/vis";

export type DeckOverlayId =
    | "grey_scatter"
    | "scatter"
    | "gate_display"
    | "selection"
    | "gate_labels";

export type ImageLayerSource = "spatial" | "ome_tiff";

export type SpatialStackEntry =
    | {
          kind: "spatial";
          layerId: string;
          visible: boolean;
          opacity?: number;
          imageSource?: ImageLayerSource;
      }
    | {
          kind: "deck";
          deckId: DeckOverlayId;
          visible: boolean;
          opacity?: number;
      };

export type SpatialLayerStackConfig = {
    stackOrder: string[];
    entries: Record<string, SpatialStackEntry>;
    spatialLayers: Record<string, LayerConfig>;
};

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

export function deckStackKey(deckId: DeckOverlayId): string {
    return `deck:${deckId}`;
}

export function isDeckStackKey(key: string): key is `deck:${DeckOverlayId}` {
    return key.startsWith("deck:");
}

export function deckIdFromStackKey(key: string): DeckOverlayId | null {
    if (!isDeckStackKey(key)) return null;
    const deckId = key.slice("deck:".length) as DeckOverlayId;
    return DECK_OVERLAY_IDS.includes(deckId) ? deckId : null;
}

function elementHasCoordinateSystem(element: unknown, coordinateSystem: string) {
    if (!element || typeof element !== "object") return false;
    if (!("coordinateSystems" in element)) return true;
    const systems = element.coordinateSystems;
    return Array.isArray(systems) && systems.includes(coordinateSystem);
}

function elementsForCoordinateSystem(
    elements: Record<string, unknown> | undefined,
    coordinateSystem: string,
) {
    return Object.entries(elements ?? {}).filter(([, element]) =>
        elementHasCoordinateSystem(element, coordinateSystem),
    );
}

export function createSpatialLayerId(type: LayerType, elementKey: string) {
    return `spatialdata-${type}-${elementKey}`;
}

export function defaultLayerConfig(
    type: LayerType,
    elementKey: string,
): LayerConfig {
    const id = createSpatialLayerId(type, elementKey);
    const base = {
        id,
        type,
        visible: true,
        opacity: 1,
        elementKey,
    };
    switch (type) {
        case "image":
            return { ...base, type: "image" };
        case "shapes":
            return {
                ...base,
                type: "shapes",
                strokeWidth: 1,
                strokeWidthUnits: "pixels",
            };
        case "points":
            return { ...base, type: "points", pointSize: 4 };
        case "labels":
            return { ...base, type: "labels" };
        default: {
            const _exhaustive: never = type;
            return _exhaustive;
        }
    }
}

export function seedSpatialLayersFromElements(
    available: ElementsByType,
    coordinateSystem: string | null,
    existingLayers: Record<string, LayerConfig> = {},
): { spatialLayers: Record<string, LayerConfig>; spatialOrder: string[] } {
    if (!coordinateSystem) {
        return { spatialLayers: { ...existingLayers }, spatialOrder: [] };
    }

    const spatialLayers = { ...existingLayers };
    const spatialOrder: string[] = [];
    const buckets: Array<{ type: LayerType; key: keyof ElementsByType }> = [
        { type: "image", key: "images" },
        { type: "shapes", key: "shapes" },
        { type: "points", key: "points" },
        { type: "labels", key: "labels" },
    ];

    for (const { type, key } of buckets) {
        const bucket = available[key];
        for (const element of bucket) {
            const id = createSpatialLayerId(type, element.key);
            if (!spatialLayers[id]) {
                spatialLayers[id] = defaultLayerConfig(type, element.key);
            }
            if (!spatialOrder.includes(id)) {
                spatialOrder.push(id);
            }
        }
    }

    return { spatialLayers, spatialOrder };
}

export function seedSpatialLayersFromSpatialData(
    spatialData: SpatialData | null,
    coordinateSystem: string | null,
    existingLayers: Record<string, LayerConfig> = {},
): { spatialLayers: Record<string, LayerConfig>; spatialOrder: string[] } {
    if (!spatialData || !coordinateSystem) {
        return { spatialLayers: { ...existingLayers }, spatialOrder: [] };
    }

    const spatialLayers = { ...existingLayers };
    const spatialOrder: string[] = [];
    const buckets: Array<{ type: LayerType; elements: Record<string, unknown> | undefined }> = [
        { type: "image", elements: spatialData.images as Record<string, unknown> },
        { type: "shapes", elements: spatialData.shapes as Record<string, unknown> },
        { type: "points", elements: spatialData.points as Record<string, unknown> },
        { type: "labels", elements: spatialData.labels as Record<string, unknown> },
    ];

    for (const { type, elements } of buckets) {
        for (const [key] of elementsForCoordinateSystem(elements, coordinateSystem)) {
            const id = createSpatialLayerId(type, key);
            if (!spatialLayers[id]) {
                spatialLayers[id] = defaultLayerConfig(type, key);
            }
            if (!spatialOrder.includes(id)) {
                spatialOrder.push(id);
            }
        }
    }

    return { spatialLayers, spatialOrder };
}

export function createDefaultSpatialLayerStack(
    spatialData: SpatialData | null,
    coordinateSystem: string | null,
): SpatialLayerStackConfig {
    const { spatialLayers, spatialOrder } = seedSpatialLayersFromSpatialData(
        spatialData,
        coordinateSystem,
    );

    const entries: Record<string, SpatialStackEntry> = {};
    const stackOrder: string[] = [];

    for (const layerId of spatialOrder) {
        entries[layerId] = {
            kind: "spatial",
            layerId,
            visible: spatialLayers[layerId]?.visible ?? true,
            opacity: spatialLayers[layerId]?.opacity ?? 1,
            imageSource: spatialLayers[layerId]?.type === "image" ? "spatial" : undefined,
        };
        stackOrder.push(layerId);
    }

    for (const deckId of DECK_OVERLAY_IDS) {
        const key = deckStackKey(deckId);
        entries[key] = {
            kind: "deck",
            deckId,
            visible: true,
            opacity: 1,
        };
        stackOrder.push(key);
    }

    return { stackOrder, entries, spatialLayers };
}

export function applySpatialLayerStack(
    target: SpatialLayerStackConfig,
    source: SpatialLayerStackConfig,
): void {
    target.stackOrder.splice(0, target.stackOrder.length, ...source.stackOrder);

    for (const key of Object.keys(target.entries)) {
        if (!(key in source.entries)) {
            delete target.entries[key];
        }
    }
    for (const [key, entry] of Object.entries(source.entries)) {
        target.entries[key] = entry;
    }

    for (const key of Object.keys(target.spatialLayers)) {
        if (!(key in source.spatialLayers)) {
            delete target.spatialLayers[key];
        }
    }
    for (const [key, layer] of Object.entries(source.spatialLayers)) {
        target.spatialLayers[key] = layer;
    }
}

export function normalizeSpatialLayerStack(
    stack: SpatialLayerStackConfig | undefined,
    spatialData: SpatialData | null,
    coordinateSystem: string | null,
): SpatialLayerStackConfig {
    if (!stack?.stackOrder?.length) {
        return createDefaultSpatialLayerStack(spatialData, coordinateSystem);
    }

    const defaults = createDefaultSpatialLayerStack(spatialData, coordinateSystem);
    const spatialLayers = { ...defaults.spatialLayers, ...stack.spatialLayers };
    const entries = { ...defaults.entries, ...stack.entries };
    const stackOrder = stack.stackOrder.filter((key) => entries[key]);

    for (const deckId of DECK_OVERLAY_IDS) {
        const key = deckStackKey(deckId);
        if (!entries[key]) {
            entries[key] = { kind: "deck", deckId, visible: true, opacity: 1 };
        }
        if (!stackOrder.includes(key)) {
            stackOrder.push(key);
        }
    }

    for (const layerId of Object.keys(spatialLayers)) {
        if (!entries[layerId]) {
            entries[layerId] = {
                kind: "spatial",
                layerId,
                visible: spatialLayers[layerId]?.visible ?? true,
                opacity: spatialLayers[layerId]?.opacity ?? 1,
            };
        }
    }

    return { stackOrder, entries, spatialLayers };
}

export function stackToCanvasState(stack: SpatialLayerStackConfig): {
    layers: Record<string, LayerConfig>;
    layerOrder: string[];
} {
    const layerOrder = stack.stackOrder.filter((key) => {
        const entry = stack.entries[key];
        return entry?.kind === "spatial" && stack.spatialLayers[entry.layerId];
    });

    const layers: Record<string, LayerConfig> = {};
    for (const key of layerOrder) {
        const entry = stack.entries[key];
        if (entry?.kind !== "spatial") continue;
        const config = stack.spatialLayers[entry.layerId];
        if (!config) continue;
        layers[config.id] = {
            ...config,
            visible: entry.visible,
            opacity: entry.opacity ?? config.opacity,
        };
    }

    return { layers, layerOrder };
}

export function canvasLayerOrder(stack: SpatialLayerStackConfig): string[] {
    return stack.stackOrder.filter((key) => {
        const entry = stack.entries[key];
        if (!entry) return false;
        if (entry.kind === "spatial") {
            return Boolean(stack.spatialLayers[entry.layerId]);
        }
        return entry.kind === "deck";
    });
}

export function entryDisplayName(
    key: string,
    stack: SpatialLayerStackConfig,
): string {
    const entry = stack.entries[key];
    if (!entry) return key;
    if (entry.kind === "deck") {
        return DECK_OVERLAY_LABELS[entry.deckId];
    }
    const layer = stack.spatialLayers[entry.layerId];
    if (!layer) return entry.layerId;
    return `${layer.type}: ${layer.elementKey}`;
}

export function isRemovableStackEntry(key: string, stack: SpatialLayerStackConfig): boolean {
    const entry = stack.entries[key];
    if (!entry) return false;
    if (entry.kind === "deck" && entry.deckId === "selection") {
        return false;
    }
    return true;
}
