import type {
    RenderStack,
    RenderStackEntry,
    RenderStackSpatialElementType,
} from "@spatialdata/layers";
import type { LayerConfig } from "@spatialdata/vis";
import { isObservable, runInAction, set } from "mobx";
import { useCallback } from "react";

import type {
    SpatialDataMdvReact,
    SpatialDataMdvReactConfig,
} from "@/react/components/SpatialDataMDVReact";
import { useChart } from "@/react/context";
import { useConfig } from "@/react/hooks";
import type { SpatialData } from "@spatialdata/core";
import {
    deckHostLayerId,
    deckIdFromHostLayerId,
    type DeckOverlayId,
} from "./host_overlay_ids";
import { spatialEntryAsLayerConfig } from "./render_stack_adapter";
import {
    createDefaultRenderStack,
    createHostOnlyRenderStack,
    insertDefaultImageLayer,
    mergeHostOverlayEntriesInPlace,
    spatialEntryId,
} from "./render_stack_defaults";
import { touchRenderStackEntry } from "./render_stack_observe";

export function reorderRenderStackEntries(
    stack: RenderStack,
    fromIndex: number,
    toIndex: number,
): void {
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
    const entry = stack.entries.find((item) => item.id === entryId);
    if (!entry) return;
    if (patch.visible !== undefined) {
        entry.visible = patch.visible;
    }
    if (patch.props) {
        for (const [key, value] of Object.entries(patch.props)) {
            if (isObservable(entry.props)) {
                set(entry.props, key, value);
            } else {
                entry.props[key] = value;
            }
        }
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
    const index = stack.entries.findIndex((entry) => entry.id === entryId);
    if (index === -1) return;
    stack.entries.splice(index, 1);
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
                ? createDefaultRenderStack(spatialData, coordinateSystem)
                : createHostOnlyRenderStack();
            chart?.finishDefaultSpatialLayerSeed();
            chart?.bumpRenderStackGeneration();
        });
        return;
    }

    runInAction(() => {
        let changed = false;
        if (chart?.seedDefaultSpatialLayers) {
            changed =
                insertDefaultImageLayer(
                    config.renderStack!,
                    spatialData,
                    coordinateSystem,
                ) || changed;
            chart.finishDefaultSpatialLayerSeed();
        }
        changed = mergeHostOverlayEntriesInPlace(config.renderStack!) || changed;
        if (changed) {
            chart?.bumpRenderStackGeneration();
        }
    });
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
 * - `chart.renderStackGeneration` bumps on commit so the canvas can refresh without
 *   replacing `config.renderStack` (which would re-render the whole dialog list).
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
    const entry = config.renderStack?.entries.find((item) => item.id === entryId);
    touchRenderStackEntry(entry);

    const patchEntry = useCallback(
        (patch: Parameters<typeof patchRenderStackEntry>[2]) => {
            runInAction(() => {
                if (!config.renderStack) return;
                patchRenderStackEntry(config.renderStack, entryId, patch);
                chart.bumpRenderStackGeneration();
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

    const layer =
        entry?.kind === "spatial" ? spatialEntryAsLayerConfig(entry) : null;

    return {
        entry,
        layer,
        patchEntry,
        patchProps,
        patchLayer,
    };
}
