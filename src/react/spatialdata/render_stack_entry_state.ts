import type { RenderStack, RenderStackEntry } from "@spatialdata/layers";
import type { LayerConfig } from "@spatialdata/vis";
import { runInAction } from "mobx";
import { useCallback } from "react";

import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "@/react/components/SpatialDataMDVReact";
import { useChart } from "@/react/context";
import { useConfig } from "@/react/hooks";
import { patchRenderStackEntry } from "./render_stack_mutations";

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

export function spatialEntryAsLayerConfig(
    entry: Extract<RenderStackEntry, { kind: "spatial" }>,
): LayerConfig {
    return {
        id: entry.id,
        type: entry.source.elementType,
        elementKey: entry.source.elementKey,
        visible: entry.visible,
        opacity: typeof entry.props.opacity === "number" ? entry.props.opacity : 1,
        ...entry.props,
    } as LayerConfig;
}

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
