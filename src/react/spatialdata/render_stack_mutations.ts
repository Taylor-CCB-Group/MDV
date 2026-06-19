import type { RenderStack, RenderStackEntry, RenderStackSpatialElementType } from "@spatialdata/layers";

import { deckHostLayerId, deckIdFromHostLayerId, type DeckOverlayId } from "./host_overlay_ids";
import { spatialEntryId } from "./render_stack_seed";

export function reorderRenderStackEntries(
    stack: RenderStack,
    fromIndex: number,
    toIndex: number,
): void {
    const next = [...stack.entries];
    const [moved] = next.splice(fromIndex, 1);
    if (!moved) return;
    next.splice(toIndex, 0, moved);
    stack.entries = next;
}

export function patchRenderStackEntry(
    stack: RenderStack,
    entryId: string,
    patch: {
        visible?: boolean;
        props?: Record<string, unknown>;
    },
): void {
    const entry = stack.entries.find((item) => item.id === entryId);
    if (!entry) return;
    if (patch.visible !== undefined) {
        entry.visible = patch.visible;
    }
    if (patch.props) {
        entry.props = { ...entry.props, ...patch.props };
    }
}

export function insertSpatialRenderStackEntry(
    stack: RenderStack,
    elementType: RenderStackSpatialElementType,
    elementKey: string,
    props: Record<string, unknown> = {},
): boolean {
    const id = spatialEntryId(elementType, elementKey);
    if (stack.entries.some((entry) => entry.id === id)) {
        return false;
    }
    stack.entries.push({
        kind: "spatial",
        id,
        visible: true,
        source: { elementType, elementKey },
        props: { opacity: 1, ...props },
    });
    return true;
}

export function insertHostRenderStackEntry(
    stack: RenderStack,
    deckId: DeckOverlayId,
): boolean {
    const id = deckHostLayerId(deckId);
    if (stack.entries.some((entry) => entry.id === id)) {
        return false;
    }
    stack.entries.push({
        kind: "host",
        id,
        visible: true,
        source: { hostLayerId: id },
        props: {},
    });
    return true;
}

export function removeRenderStackEntry(stack: RenderStack, entryId: string): void {
    stack.entries = stack.entries.filter((entry) => entry.id !== entryId);
}

export function renderStackEntryIds(stack: RenderStack): string[] {
    return stack.entries.map((entry) => entry.id);
}

export function isRemovableRenderStackEntry(entry: RenderStackEntry): boolean {
    if (entry.kind === "host") {
        return deckIdFromHostLayerId(entry.source.hostLayerId) !== null;
    }
    return entry.kind === "spatial";
}
