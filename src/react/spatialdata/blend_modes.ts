import {
    GROUP_COMPOSITING_OPERATORS,
    isGroupCompositingOperator,
    type RenderStack,
    type RenderStackBlendMode,
    type RenderStackEntry,
    type RenderStackGroupBlendMode,
    type RenderStackGroupEntry,
    resolveRenderStackGroupProps,
} from "@spatialdata/layers";
import { isObservable, remove, set } from "mobx";

import { DATA_POINTS_GROUP_ID, deckHostLayerId } from "./host_overlay_ids";

/**
 * How render-stack entries blend, as sd.js draws them (`renderStackToDeckLayers`, ADR 0007):
 *
 * - A group with any mode but pass-through, or with `isolate`, draws its layers into a target
 *   of its own and composites that onto what is beneath it, with its opacity and mode.
 * - A layer's own mode, and a group's compositing operator, act only inside such a group: they
 *   blend with what is beneath them *in that group*. Anywhere else they draw as normal.
 */

export const COLOR_BLEND_MODES = [
    "normal",
    "multiply",
    "screen",
    "additive",
    "max",
    "min",
] as const satisfies readonly RenderStackBlendMode[];

export const COMPOSITING_OPERATORS: readonly RenderStackBlendMode[] = GROUP_COMPOSITING_OPERATORS;

export const BLEND_MODE_LABELS: Record<RenderStackGroupBlendMode, string> = {
    "pass-through": "Pass through",
    normal: "Normal",
    multiply: "Multiply",
    screen: "Screen",
    additive: "Add",
    max: "Lighten (max)",
    min: "Darken (min)",
    intersect: "Intersect",
    "cut-out": "Cut out",
    xor: "Exclude (xor)",
    clip: "Clip",
};

/** What a mode does, with "this" for the layer or group it is set on. */
const BLEND_MODE_EFFECTS: Record<RenderStackGroupBlendMode, string> = {
    "pass-through": "each layer drawn as if not grouped",
    normal: "drawn over what is beneath",
    multiply: "multiplied with what is beneath, which darkens it; white leaves it unchanged",
    screen: "screened over what is beneath, which lightens it; black leaves it unchanged",
    additive: "added to what is beneath",
    max: "the lighter of this and what is beneath, per colour channel",
    min: "the darker of this and what is beneath, per colour channel",
    intersect: "what is beneath is kept only where this covers it",
    "cut-out": "what is beneath is kept only where this does not cover it",
    xor: "this and what is beneath are each kept only where the other is absent",
    clip: "this is drawn only where something is beneath it",
};

export function blendModeEffect(mode: RenderStackGroupBlendMode): string {
    return BLEND_MODE_EFFECTS[mode];
}

const GROUP_BLEND_MODES: readonly string[] = ["pass-through", ...COLOR_BLEND_MODES, ...GROUP_COMPOSITING_OPERATORS];

export function isGroupBlendMode(value: string): value is RenderStackGroupBlendMode {
    return GROUP_BLEND_MODES.includes(value);
}

export function isLeafBlendMode(value: string): value is RenderStackBlendMode {
    return value !== "pass-through" && isGroupBlendMode(value);
}

export function isCompositingOperator(mode: string): boolean {
    return isGroupCompositingOperator(mode);
}

/** The groups holding an entry, outermost first; undefined when the stack has no such entry. */
export function renderStackAncestors(stack: RenderStack, entryId: string): RenderStackGroupEntry[] | undefined {
    const path: RenderStackGroupEntry[] = [];
    const search = (entries: readonly RenderStackEntry[]): boolean => {
        for (const entry of entries) {
            if (entry.id === entryId) return true;
            if (entry.kind !== "group") continue;
            path.push(entry);
            if (search(entry.children)) return true;
            path.pop();
        }
        return false;
    };
    return search(stack.entries) ? path : undefined;
}

/**
 * Whether an entry draws inside an isolated group, at any depth: where its own blend mode, or
 * a group's compositing operator, takes effect.
 */
export function isInsideIsolatedGroup(stack: RenderStack, entryId: string): boolean {
    return (renderStackAncestors(stack, entryId) ?? []).some((group) => resolveRenderStackGroupProps(group).isolate);
}

/** A leaf entry's own blend mode, `normal` when it has none. */
export function leafBlendMode(entry: RenderStackEntry): RenderStackBlendMode {
    if (entry.kind === "group") return "normal";
    const mode = entry.props.blendMode;
    return typeof mode === "string" && isLeafBlendMode(mode) ? mode : "normal";
}

/**
 * Set a layer's own blend mode, in place. `normal` is stored as no mode, so stacks that never
 * use one carry nothing. Only the stack's structure (`useRenderStackStructure`) reads it, so
 * this needs no generation bump: the layer's config, and so its data, are untouched.
 */
export function setLeafBlendMode(stack: RenderStack, entryId: string, mode: RenderStackBlendMode): boolean {
    const entry = findLeaf(stack.entries, entryId);
    if (!entry) return false;
    if (mode === "normal") {
        if (!("blendMode" in entry.props)) return false;
        if (isObservable(entry.props)) remove(entry.props, "blendMode");
        else delete entry.props.blendMode;
        return true;
    }
    if (entry.props.blendMode === mode) return false;
    if (isObservable(entry.props)) set(entry.props, "blendMode", mode);
    else entry.props.blendMode = mode;
    return true;
}

function findLeaf(
    entries: readonly RenderStackEntry[],
    entryId: string,
): Exclude<RenderStackEntry, RenderStackGroupEntry> | undefined {
    for (const entry of entries) {
        if (entry.kind === "group") {
            const found = findLeaf(entry.children, entryId);
            if (found) return found;
        } else if (entry.id === entryId) {
            return entry;
        }
    }
    return undefined;
}

export type CompositingBadge = {
    label: string;
    title: string;
    /** Saved, but drawn as normal where the entry is (outside any isolated group). */
    inactive: boolean;
};

const NOT_ISOLATED_NOTE =
    "Drawn as Normal here: it acts only inside a group with a blend mode, or one drawn as one layer.";

/**
 * How a layer list row shows an entry's compositing, when it is not the default: a group's
 * mode or isolation, a layer's own mode, and for the data points row, their group's.
 */
export function compositingBadge(stack: RenderStack, entry: RenderStackEntry): CompositingBadge | null {
    if (entry.kind === "group") return groupBadge(stack, entry);
    if (entry.id === deckHostLayerId("scatter")) {
        const group = stack.entries.find((item) => item.id === DATA_POINTS_GROUP_ID);
        return group?.kind === "group" ? groupBadge(stack, group) : null;
    }
    const mode = leafBlendMode(entry);
    if (mode === "normal") return null;
    const inactive = !isInsideIsolatedGroup(stack, entry.id);
    return {
        label: BLEND_MODE_LABELS[mode],
        title: inactive
            ? `Blend mode: ${BLEND_MODE_LABELS[mode]}. ${NOT_ISOLATED_NOTE}`
            : `Blend mode: ${BLEND_MODE_LABELS[mode]}. In its group: ${blendModeEffect(mode)}.`,
        inactive,
    };
}

function groupBadge(stack: RenderStack, group: RenderStackGroupEntry): CompositingBadge | null {
    const { blendMode, isolate } = resolveRenderStackGroupProps(group);
    if (blendMode === "pass-through") {
        return isolate
            ? { label: "Isolated", title: "Drawn as one layer, then composited normally.", inactive: false }
            : null;
    }
    const inactive = isGroupCompositingOperator(blendMode) && !isInsideIsolatedGroup(stack, group.id);
    return {
        label: BLEND_MODE_LABELS[blendMode],
        title: inactive
            ? `Blend mode: ${BLEND_MODE_LABELS[blendMode]}. ${NOT_ISOLATED_NOTE}`
            : `Drawn as one layer, then composited: ${blendModeEffect(blendMode)}.`,
        inactive,
    };
}
