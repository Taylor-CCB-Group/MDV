import type { SpatialData, SpatialElement } from "@spatialdata/core";
import {
    RENDER_STACK_SCHEMA_VERSION,
    type RenderStack,
    type RenderStackEntry,
    type RenderStackSpatialElementType,
} from "@spatialdata/layers";

import {
    DECK_OVERLAY_IDS,
    deckHostLayerId,
} from "./host_overlay_ids";

function elementsForCoordinateSystem<T extends SpatialElement>(
    elements: Record<string, T> | undefined,
    coordinateSystem: string,
): Array<[string, T]> {
    return Object.entries(elements ?? {}).filter(([, element]) =>
        element.coordinateSystems.includes(coordinateSystem),
    );
}

export function spatialEntryId(
    elementType: RenderStackSpatialElementType,
    elementKey: string,
) {
    const id = Math.random().toString(36).substring(2, 4+2);
    return `spatialdata-${elementType}-${elementKey}-#${id}`;
}

function spatialEntry(
    elementType: RenderStackSpatialElementType,
    elementKey: string,
    props: Record<string, unknown> = {},
): RenderStackEntry {
    return {
        kind: "spatial",
        id: spatialEntryId(elementType, elementKey),
        visible: true,
        source: { elementType, elementKey },
        props: { opacity: 1, ...props },
    };
}

function hostEntry(deckId: (typeof DECK_OVERLAY_IDS)[number]): RenderStackEntry {
    const hostLayerId = deckHostLayerId(deckId);
    return {
        kind: "host",
        id: hostLayerId,
        visible: true,
        source: { hostLayerId },
        props: {},
    };
}

export function defaultPropsForSpatialElement(
    elementType: RenderStackSpatialElementType,
): Record<string, unknown> {
    switch (elementType) {
        case "shapes":
            return {
                fillColor: [100, 149, 237, 180],
                strokeWidth: 1,
                strokeWidthUnits: "pixels",
            };
        case "points":
            return { pointSize: 4 };
        default:
            return {};
    }
}

function addAvailableSpatialEntries<T extends SpatialElement>(
    spatialEntries: RenderStackEntry[],
    elementType: RenderStackSpatialElementType,
    elements: Record<string, T> | undefined,
    coordinateSystem: string,
) {
    for (const [key] of elementsForCoordinateSystem(elements, coordinateSystem)) {
        spatialEntries.push(
            spatialEntry(elementType, key, defaultPropsForSpatialElement(elementType)),
        );
    }
}

/** All spatial elements in the store — for dialog "insert layer" options only. */
export function listAvailableSpatialEntries(
    spatialData: SpatialData,
    coordinateSystem: string,
): RenderStackEntry[] {
    const spatialEntries: RenderStackEntry[] = [];
    addAvailableSpatialEntries(spatialEntries, "image", spatialData.images, coordinateSystem);
    addAvailableSpatialEntries(spatialEntries, "shapes", spatialData.shapes, coordinateSystem);
    addAvailableSpatialEntries(spatialEntries, "points", spatialData.points, coordinateSystem);
    addAvailableSpatialEntries(spatialEntries, "labels", spatialData.labels, coordinateSystem);

    return spatialEntries;
}

function firstImageSpatialEntry(
    spatialData: SpatialData,
    coordinateSystem: string,
): RenderStackEntry | null {
    const matching = elementsForCoordinateSystem(spatialData.images, coordinateSystem);
    const firstKey = matching[0]?.[0];
    if (!firstKey) return null;
    return spatialEntry("image", firstKey);
}

export function createHostOnlyRenderStack(): RenderStack {
    return {
        schemaVersion: RENDER_STACK_SCHEMA_VERSION,
        entries: DECK_OVERLAY_IDS.map((deckId) => hostEntry(deckId)),
    };
}

export function createDefaultRenderStack(
    spatialData: SpatialData,
    coordinateSystem: string,
): RenderStack {
    const spatialEntries: RenderStackEntry[] = [];
    const imageEntry = firstImageSpatialEntry(spatialData, coordinateSystem);
    if (imageEntry) {
        spatialEntries.push(imageEntry);
    }
    const hostEntries = DECK_OVERLAY_IDS.map((deckId) => hostEntry(deckId));
    return {
        schemaVersion: RENDER_STACK_SCHEMA_VERSION,
        entries: [...spatialEntries, ...hostEntries],
    };
}

/** Insert the first available image layer when the stack has no image entries yet. */
export function insertDefaultImageLayer(
    stack: RenderStack,
    spatialData: SpatialData,
    coordinateSystem: string,
): boolean {
    const hasImage = stack.entries.some(
        (entry) => entry.kind === "spatial" && entry.source.elementType === "image",
    );
    if (hasImage) return false;

    const imageEntry = firstImageSpatialEntry(spatialData, coordinateSystem);
    if (!imageEntry) return false;

    const firstHostIndex = stack.entries.findIndex((entry) => entry.kind === "host");
    if (firstHostIndex === -1) {
        stack.entries.push(imageEntry);
    } else {
        stack.entries.splice(firstHostIndex, 0, imageEntry);
    }
    return true;
}

/**
 * Ensure required host overlay entries exist. Never injects spatial layers —
 * persisted stacks with zero spatial layers stay as saved.
 */
export function normalizeRenderStack(stack: RenderStack | undefined): RenderStack {
    if (!stack?.entries?.length) {
        return createHostOnlyRenderStack();
    }

    const entriesById = new Map<string, RenderStackEntry>();
    for (const entry of stack.entries) {
        entriesById.set(entry.id, entry);
    }
    for (const deckId of DECK_OVERLAY_IDS) {
        const id = deckHostLayerId(deckId);
        if (!entriesById.has(id)) {
            entriesById.set(id, hostEntry(deckId));
        }
    }

    const order = stack.entries.map((entry) => entry.id);
    for (const deckId of DECK_OVERLAY_IDS) {
        const id = deckHostLayerId(deckId);
        if (!order.includes(id)) {
            order.push(id);
        }
    }

    return {
        schemaVersion: RENDER_STACK_SCHEMA_VERSION,
        entries: order
            .map((id) => entriesById.get(id))
            .filter((entry): entry is RenderStackEntry => entry !== undefined),
    };
}

export function mergeHostOverlayEntriesInPlace(stack: RenderStack): boolean {
    const normalized = normalizeRenderStack(stack);
    const beforeIds = stack.entries.map((entry) => entry.id).join("\0");
    const afterIds = normalized.entries.map((entry) => entry.id).join("\0");
    if (beforeIds === afterIds) {
        return false;
    }

    const existingById = new Map(stack.entries.map((entry) => [entry.id, entry]));
    const mergedEntries = normalized.entries.map(
        (entry) => existingById.get(entry.id) ?? entry,
    );
    stack.schemaVersion = normalized.schemaVersion;
    stack.entries.splice(0, stack.entries.length, ...mergedEntries);
    return true;
}
