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
    return `spatialdata-${elementType}-${elementKey}`;
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

export function normalizeRenderStack(
    stack: RenderStack | undefined,
    spatialData: SpatialData,
    coordinateSystem: string,
): RenderStack {
    const defaults = createDefaultRenderStack(spatialData, coordinateSystem);
    if (!stack?.entries?.length) {
        return defaults;
    }

    const entriesById = new Map<string, RenderStackEntry>();
    for (const entry of defaults.entries) {
        entriesById.set(entry.id, entry);
    }
    for (const entry of stack.entries) {
        entriesById.set(entry.id, entry);
    }

    const order = stack.entries.map((entry) => entry.id);
    for (const entry of defaults.entries) {
        if (!order.includes(entry.id)) {
            order.push(entry.id);
        }
    }

    return {
        schemaVersion: RENDER_STACK_SCHEMA_VERSION,
        entries: order
            .map((id) => entriesById.get(id))
            .filter((entry): entry is RenderStackEntry => entry !== undefined),
    };
}
