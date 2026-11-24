import type { LayerChannelConfig } from "@spatialdata/avivatorish";
import type { SpatialData } from "@spatialdata/core";
import {
    findRenderStackEntry,
    type LabelsChannelProps,
    type RenderStack,
    type RenderStackEntry,
    type RenderStackGroupBlendMode,
    type RenderStackGroupEntry,
    type RenderStackSpatialElementType,
} from "@spatialdata/layers";
import type { LayerConfig } from "@spatialdata/vis";
import { isObservable, runInAction, set } from "mobx";
import { useCallback } from "react";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "@/react/components/SpatialDataMDVReact";
import { useChart } from "@/react/context";
import { useConfig } from "@/react/hooks";
import { getRandomString } from "@/utilities/Utilities";
import {
    DATA_POINTS_GROUP_ID,
    DATA_POINTS_GROUP_MEMBERS,
    type DeckOverlayId,
    deckHostLayerId,
    deckIdFromHostLayerId,
} from "./host_overlay_ids";

export { DATA_POINTS_GROUP_ID };

import { spatialEntryAsLayerConfig } from "./render_stack_adapter";
import {
    createDefaultRenderStack,
    createHostOnlyRenderStack,
    insertDefaultImageLayer,
    isOverlayEntry,
    mergeHostOverlayEntriesInPlace,
    spatialEntryId,
} from "./render_stack_defaults";
import { touchRenderStackEntry } from "./render_stack_observe";

/** The image channel fields, which `useLayerChannelState` writes. */
const IMAGE_CHANNEL_CONFIG_KEYS = [
    "channelIds",
    "colors",
    "contrastLimits",
    "channelsVisible",
    "selections",
] as const satisfies readonly (keyof LayerChannelConfig)[];

/** A labels layer's per-channel appearance, which shares the `channels` prop with image layers. */
const LABELS_CHANNEL_CONFIG_KEYS = [
    "channelOpacities",
    "channelOutlineOpacities",
    "channelsFilled",
    "channelStrokeWidths",
] as const satisfies readonly (keyof LabelsChannelProps)[];

const CHANNEL_CONFIG_KEYS = [...IMAGE_CHANNEL_CONFIG_KEYS, ...LABELS_CHANNEL_CONFIG_KEYS];

type ChannelConfigPatch = LayerChannelConfig &
    Partial<Pick<LabelsChannelProps, (typeof LABELS_CHANNEL_CONFIG_KEYS)[number]>>;

function patchRecordInPlace(existing: Record<string, unknown>, next: Record<string, unknown>): void {
    for (const [key, value] of Object.entries(next)) {
        if (isObservable(existing)) {
            set(existing, key, value);
        } else {
            existing[key] = value;
        }
    }
}

function patchChannelsInPlace(existing: Record<string, unknown>, next: ChannelConfigPatch): void {
    for (const key of CHANNEL_CONFIG_KEYS) {
        const value = next[key];
        if (value === undefined) continue;
        if (isObservable(existing)) {
            set(existing, key, value);
        } else {
            existing[key] = value;
        }
    }
}

export function patchRenderStackEntryPropsInPlace(
    props: Record<string, unknown>,
    patch: Record<string, unknown>,
): void {
    for (const [key, value] of Object.entries(patch)) {
        if (
            key === "channels" &&
            value &&
            typeof value === "object" &&
            props.channels &&
            typeof props.channels === "object"
        ) {
            patchChannelsInPlace(props.channels as Record<string, unknown>, value as ChannelConfigPatch);
            continue;
        }
        if (
            key === "vivLayerProps" &&
            value &&
            typeof value === "object" &&
            props.vivLayerProps &&
            typeof props.vivLayerProps === "object"
        ) {
            patchRecordInPlace(props.vivLayerProps as Record<string, unknown>, value as Record<string, unknown>);
            continue;
        }
        if (isObservable(props)) {
            set(props, key, value);
        } else {
            props[key] = value;
        }
    }
}

export function reorderRenderStackEntries(stack: RenderStack, fromIndex: number, toIndex: number): void {
    const [moved] = stack.entries.splice(fromIndex, 1);
    if (!moved) return;
    stack.entries.splice(toIndex, 0, moved);
}

export function patchRenderStackEntry(
    stack: RenderStack,
    entryId: string,
    patch: {
        visible?: boolean;
        props?: Record<string, unknown>;
    },
): void {
    const entry = findRenderStackEntry(stack, entryId);
    if (!entry) return;
    if (patch.visible !== undefined) {
        entry.visible = patch.visible;
    }
    if (patch.props) {
        patchRenderStackEntryPropsInPlace(entry.props, patch.props);
    }
}

/**
 * Add a new entry on top of `parentId`'s children, or, for null or a group that is not
 * there, on top of the top-level entries below the MDV overlays: the overlays always
 * draw above SpatialData layers.
 */
function insertNewEntry(stack: RenderStack, entry: RenderStackEntry, parentId: string | null): void {
    const parent = parentId === null ? undefined : findRenderStackEntry(stack, parentId);
    if (parent?.kind === "group") {
        parent.children.push(entry);
        return;
    }
    const firstOverlayIndex = stack.entries.findIndex(isOverlayEntry);
    stack.entries.splice(firstOverlayIndex === -1 ? stack.entries.length : firstOverlayIndex, 0, entry);
}

/** Add a layer for a SpatialData element (see `insertNewEntry`). Returns its id. */
export function insertSpatialRenderStackEntry(
    stack: RenderStack,
    elementType: RenderStackSpatialElementType,
    elementKey: string,
    props: Record<string, unknown> = {},
    parentId: string | null = null,
): string | null {
    const id = spatialEntryId(elementType, elementKey);
    if (findRenderStackEntry(stack, id)) {
        return null;
    }
    insertNewEntry(
        stack,
        { kind: "spatial", id, visible: true, source: { elementType, elementKey }, props: { opacity: 1, ...props } },
        parentId,
    );
    return id;
}

export function insertHostRenderStackEntry(stack: RenderStack, deckId: DeckOverlayId): boolean {
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

/** Where an entry sits: the list holding it (the stack's top level or a group's children). */
export function findRenderStackSiblings(
    stack: RenderStack,
    entryId: string,
): { siblings: RenderStackEntry[]; index: number; parent: RenderStackGroupEntry | null } | undefined {
    const search = (
        siblings: RenderStackEntry[],
        parent: RenderStackGroupEntry | null,
    ): ReturnType<typeof findRenderStackSiblings> => {
        const index = siblings.findIndex((entry) => entry.id === entryId);
        if (index !== -1) return { siblings, index, parent };
        for (const entry of siblings) {
            if (entry.kind !== "group") continue;
            const found = search(entry.children, entry);
            if (found) return found;
        }
        return undefined;
    };
    return search(stack.entries, null);
}

export function removeRenderStackEntry(stack: RenderStack, entryId: string): void {
    const found = findRenderStackSiblings(stack, entryId);
    found?.siblings.splice(found.index, 1);
}

export function renderStackGroupId(): string {
    return `group-#${getRandomString()}`;
}

/** Add an empty group where a new layer goes (see `insertNewEntry`). Returns its id. */
export function insertRenderStackGroup(stack: RenderStack, parentId: string | null = null): string {
    const id = renderStackGroupId();
    insertNewEntry(stack, { kind: "group", id, visible: true, props: {}, children: [] }, parentId);
    return id;
}

/** Replace a group with its children, in its place. */
export function ungroupRenderStackEntry(stack: RenderStack, groupId: string): boolean {
    const found = findRenderStackSiblings(stack, groupId);
    const group = found?.siblings[found.index];
    if (!found || group?.kind !== "group") return false;
    found.siblings.splice(found.index, 1, ...group.children);
    return true;
}

/**
 * Move an entry into `parentId` (a group, or the top level for null), just after
 * `afterId` there, or first when `afterId` is null. Refuses to put a group inside itself.
 */
export function moveRenderStackEntry(
    stack: RenderStack,
    entryId: string,
    parentId: string | null,
    afterId: string | null,
): boolean {
    const from = findRenderStackSiblings(stack, entryId);
    if (!from) return false;
    const moved = from.siblings[from.index];
    let target = stack.entries;
    if (parentId !== null) {
        const parent = findRenderStackEntry(stack, parentId);
        if (parent?.kind !== "group") return false;
        if (moved.kind === "group" && findRenderStackEntry({ ...stack, entries: [moved] }, parentId)) return false;
        target = parent.children;
    }
    from.siblings.splice(from.index, 1);
    const afterIndex = afterId === null ? -1 : target.findIndex((entry) => entry.id === afterId);
    target.splice(afterIndex + 1, 0, moved);
    if (target === stack.entries && !isOverlayEntry(moved)) {
        // The dialog lists only spatial rows, so `afterId` can sit below an overlay that a
        // stack saved with them interleaved still has among them. Keep every overlay above
        // the spatial layers, in its own order.
        const overlays = stack.entries.filter(isOverlayEntry);
        const spatial = stack.entries.filter((entry) => !isOverlayEntry(entry));
        stack.entries.splice(0, stack.entries.length, ...spatial, ...overlays);
    }
    return true;
}

export function renderStackEntryIds(stack: RenderStack): string[] {
    return stack.entries.map((entry) => entry.id);
}

/**
 * Ensure host overlays exist and, only for brand-new charts, insert the default image layer.
 * Persisted stacks (including host-only with zero spatial layers) are left unchanged.
 */
export function seedRenderStackFromSpatialData(
    config: SpatialDataMdvReactConfig,
    spatialData: SpatialData,
    coordinateSystem: string,
    chart?: SpatialDataMdvReact,
): void {
    if (!config.renderStack) {
        runInAction(() => {
            config.renderStack = chart?.seedDefaultSpatialLayers
                ? createDefaultRenderStack(spatialData, coordinateSystem, chart.defaultImageLayerProps)
                : createHostOnlyRenderStack();
            chart?.finishDefaultSpatialLayerSeed();
            chart?.bumpRenderStackGeneration();
        });
        return;
    }

    runInAction(() => {
        let changed = false;
        if (!config.renderStack) {
            throw new Error(`Unexpected config with no renderStack ${config}`);
        }
        if (chart?.seedDefaultSpatialLayers) {
            changed =
                insertDefaultImageLayer(
                    config.renderStack,
                    spatialData,
                    coordinateSystem,
                    chart.defaultImageLayerProps,
                ) || changed;
            chart.finishDefaultSpatialLayerSeed();
        }
        changed = mergeHostOverlayEntriesInPlace(config.renderStack) || changed;
        if (changed) {
            chart?.bumpRenderStackGeneration();
        } else {
            chart?.bumpRenderStackPropsGeneration();
        }
    });
}

/**
 * The group the MDV data points (the grey and coloured scatter) are composited in. It
 * exists only while their compositing differs from the default, so stacks that never
 * use it carry no group and pay nothing for one.
 */

export type CompositingProps = {
    opacity: number;
    blendMode: RenderStackGroupBlendMode;
    isolate: boolean;
};

export const DEFAULT_COMPOSITING: CompositingProps = { opacity: 1, blendMode: "pass-through", isolate: false };

export function isDefaultCompositing(props: CompositingProps): boolean {
    return (
        props.opacity === DEFAULT_COMPOSITING.opacity &&
        props.blendMode === DEFAULT_COMPOSITING.blendMode &&
        props.isolate === DEFAULT_COMPOSITING.isolate
    );
}

export function getDataPointsGroup(stack: RenderStack): RenderStackGroupEntry | undefined {
    const entry = stack.entries.find((item) => item.id === DATA_POINTS_GROUP_ID);
    return entry?.kind === "group" ? entry : undefined;
}

/**
 * Set how the data points are composited. Default compositing takes the scatter layers
 * back out of their group; anything else puts them in it, at the lower one's place.
 * Returns whether the stack's structure changed; editing an existing group's props
 * patches it in place and returns false.
 */
export function setDataPointsCompositing(stack: RenderStack, props: CompositingProps): boolean {
    const group = getDataPointsGroup(stack);
    if (isDefaultCompositing(props)) {
        if (!group) return false;
        stack.entries.splice(stack.entries.indexOf(group), 1, ...group.children);
        return true;
    }
    if (group) {
        Object.assign(group.props, props);
        return false;
    }
    const members = stack.entries.filter((entry) => DATA_POINTS_GROUP_MEMBERS.includes(entry.id));
    if (members.length === 0) return false;
    const at = stack.entries.indexOf(members[0]);
    for (const member of members) stack.entries.splice(stack.entries.indexOf(member), 1);
    stack.entries.splice(at, 0, {
        kind: "group",
        id: DATA_POINTS_GROUP_ID,
        visible: true,
        props: { ...props },
        children: members,
    });
    return true;
}

export function isRemovableRenderStackEntry(entry: RenderStackEntry): boolean {
    if (entry.kind === "host") {
        return deckIdFromHostLayerId(entry.source.hostLayerId) !== null;
    }
    return entry.kind === "spatial";
}

/**
 * MobX render-stack control model for the layer dialog and viewer.
 *
 * - `config.renderStack` stays a stable object; entries and props are patched in place.
 * - `chart.renderStackGeneration` bumps only for entry-level changes (e.g. visibility)
 *   that need a coarse canvas refresh. In-place `props` edits bump
 *   `renderStackPropsGeneration` so cosmetic channel/tone/opacity/annotation edits
 *   refresh viewer inputs without resetting spatial renderer passthrough.
 * - UI rows use `useRenderStackEntry(entryId)` inside `observer` components so each row
 *   subscribes only to its own entry. Panels receive plain derived props at the boundary.
 *
 * Call hooks only from `observer` components.
 */
export function useRenderStackMutation() {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();

    return useCallback(
        (mutator: (stack: RenderStack) => void) => {
            runInAction(() => {
                if (!config.renderStack) return;
                mutator(config.renderStack);
                chart.bumpRenderStackGeneration();
            });
        },
        [chart, config],
    );
}

export function useRenderStackEntry(entryId: string) {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const entry = config.renderStack && findRenderStackEntry(config.renderStack, entryId);
    touchRenderStackEntry(entry);

    const patchEntry = useCallback(
        (patch: Parameters<typeof patchRenderStackEntry>[2]) => {
            runInAction(() => {
                if (!config.renderStack) return;
                patchRenderStackEntry(config.renderStack, entryId, patch);
                // Props-only in-place edits are observed via touchRenderStackEntry /
                // renderStackSpatialRevision plus a lightweight prop token; avoid the
                // coarse generation bump that recreates viv passthrough state.
                if (patch.props) {
                    chart.bumpRenderStackPropsGeneration();
                }
                if (patch.visible !== undefined) {
                    chart.bumpRenderStackGeneration();
                }
            });
        },
        [chart, config, entryId],
    );

    const patchProps = useCallback(
        (props: Record<string, unknown>) => {
            patchEntry({ props });
        },
        [patchEntry],
    );

    const patchLayer = useCallback(
        (updates: Partial<LayerConfig>) => {
            const { visible, ...rest } = updates as Partial<LayerConfig> & {
                visible?: boolean;
            };
            const patch: Parameters<typeof patchRenderStackEntry>[2] = {};
            if (visible !== undefined) patch.visible = visible;
            if (Object.keys(rest).length > 0) {
                patch.props = rest as Record<string, unknown>;
            }
            patchEntry(patch);
        },
        [patchEntry],
    );

    const layer = entry?.kind === "spatial" ? spatialEntryAsLayerConfig(entry) : null;

    return {
        entry,
        layer,
        patchEntry,
        patchProps,
        patchLayer,
    };
}
