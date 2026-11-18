import type { LayerChannelConfig } from "@spatialdata/avivatorish";
import type { SpatialData, SpatialElement } from "@spatialdata/core";
import {
    findRenderStackEntry,
    migrateRenderStack,
    RENDER_STACK_SCHEMA_VERSION,
    type RenderStack,
    type RenderStackEntry,
    type RenderStackSpatialElementType,
    type RenderStackV1,
    walkRenderStackEntries,
} from "@spatialdata/layers";

import { getRandomString } from "@/utilities/Utilities";

import { DECK_OVERLAY_IDS, DECK_OVERLAY_VISIBILITY_FOLLOWERS, deckHostLayerId } from "./host_overlay_ids";

function elementsForCoordinateSystem<T extends SpatialElement>(
    elements: Record<string, T> | undefined,
    coordinateSystem: string,
): Array<[string, T]> {
    return Object.entries(elements ?? {}).filter(([, element]) => element.coordinateSystems.includes(coordinateSystem));
}

export function spatialEntryId(elementType: RenderStackSpatialElementType, elementKey: string) {
    // The same element can appear in multiple stack entries (different roles /
    // blend modes), so each entry needs its own id. The random suffix only has to
    // be unique within a single stack — a smaller scope than the project-wide
    // chart ids that already use the 6-char `getRandomString()` — so that length
    // is ample here while keeping the id readable in saved JSON.
    return `spatialdata-${elementType}-${elementKey}-#${getRandomString()}`;
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

export function defaultPropsForSpatialElement(elementType: RenderStackSpatialElementType): Record<string, unknown> {
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
        spatialEntries.push(spatialEntry(elementType, key, defaultPropsForSpatialElement(elementType)));
    }
}

/** All spatial elements in the store — for dialog "insert layer" options only. */
export function listAvailableSpatialEntries(spatialData: SpatialData, coordinateSystem: string): RenderStackEntry[] {
    const spatialEntries: RenderStackEntry[] = [];
    addAvailableSpatialEntries(spatialEntries, "image", spatialData.images, coordinateSystem);
    addAvailableSpatialEntries(spatialEntries, "shapes", spatialData.shapes, coordinateSystem);
    addAvailableSpatialEntries(spatialEntries, "points", spatialData.points, coordinateSystem);
    addAvailableSpatialEntries(spatialEntries, "labels", spatialData.labels, coordinateSystem);

    return spatialEntries;
}

/**
 * Props for the image layer a chart adds when it starts without one, such as the channel settings
 * of a Viv chart converted to the SpatialData.js viewer.
 */
export type DefaultImageLayerProps = {
    channels?: LayerChannelConfig;
    vivLayerProps?: Record<string, unknown>;
};

function firstImageSpatialEntry(
    spatialData: SpatialData,
    coordinateSystem: string,
    props?: DefaultImageLayerProps,
): RenderStackEntry | null {
    const matching = elementsForCoordinateSystem(spatialData.images, coordinateSystem);
    const firstKey = matching[0]?.[0];
    if (!firstKey) return null;
    // a copy, so that editing the layer does not edit `props`
    return spatialEntry("image", firstKey, props ? structuredClone(props) : {});
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
    imageProps?: DefaultImageLayerProps,
): RenderStack {
    const spatialEntries: RenderStackEntry[] = [];
    const imageEntry = firstImageSpatialEntry(spatialData, coordinateSystem, imageProps);
    if (imageEntry) {
        spatialEntries.push(imageEntry);
    }
    const hostEntries = DECK_OVERLAY_IDS.map((deckId) => hostEntry(deckId));
    return {
        schemaVersion: RENDER_STACK_SCHEMA_VERSION,
        entries: [...spatialEntries, ...hostEntries],
    };
}

/**
 * An MDV overlay, or a group holding them (the data points' group). Overlays always draw
 * above SpatialData layers, so new spatial layers go below the first of these.
 */
export function isOverlayEntry(entry: RenderStackEntry): boolean {
    if (entry.kind === "host") return true;
    return entry.kind === "group" && Array.from(walkRenderStackEntries(entry.children)).some((e) => e.kind === "host");
}

/** Insert the first available image layer when the stack has no image entries yet. */
export function insertDefaultImageLayer(
    stack: RenderStack,
    spatialData: SpatialData,
    coordinateSystem: string,
    imageProps?: DefaultImageLayerProps,
): boolean {
    const hasImage = Array.from(walkRenderStackEntries(stack.entries)).some(
        (entry) => entry.kind === "spatial" && entry.source.elementType === "image",
    );
    if (hasImage) return false;

    const imageEntry = firstImageSpatialEntry(spatialData, coordinateSystem, imageProps);
    if (!imageEntry) return false;

    const firstOverlayIndex = stack.entries.findIndex(isOverlayEntry);
    stack.entries.splice(firstOverlayIndex === -1 ? stack.entries.length : firstOverlayIndex, 0, imageEntry);
    return true;
}

/**
 * A saved stack as the current schema. Saved projects may hold a v1 stack, which
 * `RenderStack` can't describe, so this is the first thing to touch one on load.
 * Notices are child ids a v1 group named that the migration had to drop.
 */
export function migrateSavedRenderStack(saved: RenderStackV1 | RenderStack): RenderStack {
    const { stack, notices } = migrateRenderStack(saved);
    for (const notice of notices) {
        console.warn(
            `Render stack migration dropped '${notice.childId}' from group '${notice.groupId}' (${notice.reason})`,
        );
    }
    for (const id of dropDuplicateEntries(stack)) {
        console.warn(`Render stack had '${id}' more than once; kept the first`);
    }
    return stack;
}

/**
 * Remove every entry whose id an earlier one (depth-first, in draw order) already has,
 * in place, and return their ids.
 *
 * Builds from before groups (sd.js 0.15 and older) read only the top level, and add any
 * MDV overlay they don't find there. A stack they save after the data points were put
 * in their group therefore holds those overlays twice: in the group, and again at the
 * top level, after it.
 */
function dropDuplicateEntries(stack: RenderStack): string[] {
    const seen = new Set<string>();
    const dropped: string[] = [];
    const prune = (list: RenderStackEntry[]) => {
        for (let i = 0; i < list.length; ) {
            const entry = list[i];
            if (seen.has(entry.id)) {
                dropped.push(entry.id);
                list.splice(i, 1);
                continue;
            }
            seen.add(entry.id);
            if (entry.kind === "group") prune(entry.children);
            i++;
        }
    };
    prune(stack.entries);
    return dropped;
}

/** `migrateSavedRenderStack` for a stack that is already observable: patched in place. */
function migrateRenderStackInPlace(stack: RenderStack): boolean {
    const countBefore = Array.from(walkRenderStackEntries(stack.entries)).length;
    const migrated = migrateSavedRenderStack(stack);
    if (migrated === stack) return Array.from(walkRenderStackEntries(stack.entries)).length !== countBefore;
    stack.schemaVersion = migrated.schemaVersion;
    stack.entries.splice(0, stack.entries.length, ...migrated.entries);
    return true;
}

/**
 * Ensure required host overlay entries exist, at any depth. A missing one goes at the
 * top level, just above the entry that holds the overlay it follows in drawing order.
 * Never injects spatial layers — persisted stacks with zero spatial layers stay as saved.
 */
export function normalizeRenderStack(stack: RenderStack | undefined): RenderStack {
    if (!stack?.entries?.length) {
        return createHostOnlyRenderStack();
    }
    const migrated = migrateSavedRenderStack(stack);

    const entriesById = new Map<string, RenderStackEntry>();
    // The top-level entry each entry sits under (itself, at the top level).
    const topLevelOf = new Map<string, string>();
    for (const top of migrated.entries) {
        for (const entry of walkRenderStackEntries([top])) {
            entriesById.set(entry.id, entry);
            topLevelOf.set(entry.id, top.id);
        }
    }
    const added = new Set<string>();
    for (const deckId of DECK_OVERLAY_IDS) {
        const id = deckHostLayerId(deckId);
        if (!entriesById.has(id)) {
            entriesById.set(id, hostEntry(deckId));
            added.add(id);
        }
    }
    // An added follower (e.g. the grey scatter in a stack saved before it existed) takes its
    // leader's visibility, since the layer dialog offers no checkbox of its own for it.
    for (const leaderId of DECK_OVERLAY_IDS) {
        const leader = entriesById.get(deckHostLayerId(leaderId));
        if (!leader) continue;
        for (const followerId of DECK_OVERLAY_VISIBILITY_FOLLOWERS[leaderId] ?? []) {
            const follower = entriesById.get(deckHostLayerId(followerId));
            if (follower && added.has(follower.id)) follower.visible = leader.visible;
        }
    }

    // A missing overlay goes just above the overlay it follows in drawing order, so one
    // added to MDV later (e.g. density, above the data points) is not drawn over the rest.
    // An overlay inside a group stands for the group's top-level position.
    const order = migrated.entries.map((entry) => entry.id);
    const positionOf = (deckId: (typeof DECK_OVERLAY_IDS)[number]) => {
        const id = deckHostLayerId(deckId);
        return order.indexOf(topLevelOf.get(id) ?? id);
    };
    DECK_OVERLAY_IDS.forEach((deckId, index) => {
        const id = deckHostLayerId(deckId);
        if (!added.has(id)) return;
        const below = DECK_OVERLAY_IDS.slice(0, index)
            .map(positionOf)
            .filter((position) => position >= 0);
        if (below.length > 0) {
            order.splice(Math.max(...below) + 1, 0, id);
            return;
        }
        const above = DECK_OVERLAY_IDS.slice(index + 1)
            .map(positionOf)
            .filter((position) => position >= 0);
        if (above.length > 0) order.splice(Math.min(...above), 0, id);
        else order.push(id);
    });

    return {
        schemaVersion: RENDER_STACK_SCHEMA_VERSION,
        entries: order
            .map((id) => entriesById.get(id))
            .filter((entry): entry is RenderStackEntry => entry !== undefined),
    };
}

export function mergeHostOverlayEntriesInPlace(stack: RenderStack): boolean {
    const migrated = migrateRenderStackInPlace(stack);
    const normalized = normalizeRenderStack(stack);
    const beforeIds = stack.entries.map((entry) => entry.id).join("\0");
    const afterIds = normalized.entries.map((entry) => entry.id).join("\0");
    if (beforeIds === afterIds) {
        return migrated;
    }

    const existingById = new Map(stack.entries.map((entry) => [entry.id, entry]));
    const mergedEntries = normalized.entries.map((entry) => existingById.get(entry.id) ?? entry);
    stack.schemaVersion = normalized.schemaVersion;
    stack.entries.splice(0, stack.entries.length, ...mergedEntries);
    return true;
}

/**
 * Give each following overlay its leader's visibility (the grey scatter follows the data
 * points). The layer dialog only offers the leader's checkbox, so a saved stack where the
 * two disagree would have no way to fix it. For a stack being restored only: during a
 * session a follower can be turned off on its own by a render failure, which this would undo.
 */
export function syncOverlayFollowerVisibility(stack: RenderStack): void {
    for (const leaderId of DECK_OVERLAY_IDS) {
        const followers = DECK_OVERLAY_VISIBILITY_FOLLOWERS[leaderId];
        const leader = findRenderStackEntry(stack, deckHostLayerId(leaderId));
        if (!followers || !leader) continue;
        for (const followerId of followers) {
            const follower = findRenderStackEntry(stack, deckHostLayerId(followerId));
            if (follower) follower.visible = leader.visible;
        }
    }
}
