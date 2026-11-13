import type { RenderStack, RenderStackEntry } from "@spatialdata/layers";

import { DECK_OVERLAY_LABELS, deckIdFromHostLayerId } from "./host_overlay_ids";

/**
 * Where a user-chosen layer name sits on an entry's `props`.
 *
 * `props` is the only open part of a render-stack entry — its top level is a strict
 * schema — so the name lives there, as MDV's own field. It means nothing to the viewer
 * and is taken off at the viewer boundary (`withoutMdvOnlyProps`).
 */
export const DISPLAY_NAME_PROP = "displayName";

/** The name the user gave this entry, if any; blank names count as none. */
export function customDisplayName(entry: RenderStackEntry): string | undefined {
    const value = entry.props?.[DISPLAY_NAME_PROP];
    if (typeof value !== "string") return undefined;
    const trimmed = value.trim();
    return trimmed === "" ? undefined : trimmed;
}

/** What an entry is called before anyone names it: what it draws. */
export function defaultRenderStackEntryName(entry: RenderStackEntry): string {
    if (entry.kind === "host") {
        const deckId = deckIdFromHostLayerId(entry.source.hostLayerId);
        if (deckId) {
            return DECK_OVERLAY_LABELS[deckId];
        }
        return entry.source.hostLayerId;
    }
    if (entry.kind === "spatial") {
        return `${entry.source.elementType}: ${entry.source.elementKey}`;
    }
    return entry.id;
}

function sameSpatialElement(a: RenderStackEntry, b: RenderStackEntry): boolean {
    return (
        a.kind === "spatial" &&
        b.kind === "spatial" &&
        a.source.elementType === b.source.elementType &&
        a.source.elementKey === b.source.elementKey
    );
}

/**
 * The name to show for an entry.
 *
 * A user-chosen name wins. Otherwise, with `stack` given, an unnamed entry that shares
 * its element with other entries is numbered in stack order (`shapes: cells (2)`), so
 * two layers drawing the same element can be told apart before anyone renames them.
 */
export function renderStackEntryDisplayName(entry: RenderStackEntry, stack?: RenderStack): string {
    return customDisplayName(entry) ?? numberedDefaultName(entry, stack);
}

/** The default name, numbered among entries of the same element (see `renderStackEntryDisplayName`). */
export function numberedDefaultName(entry: RenderStackEntry, stack?: RenderStack): string {
    const name = defaultRenderStackEntryName(entry);
    if (!stack || entry.kind !== "spatial") return name;
    const siblings = stack.entries.filter((item) => sameSpatialElement(item, entry));
    if (siblings.length < 2) return name;
    return `${name} (${siblings.indexOf(entry) + 1})`;
}

export function renderStackOrderLabel(stack: RenderStack): string {
    return stack.entries.map((entry) => renderStackEntryDisplayName(entry, stack)).join(" → ");
}
