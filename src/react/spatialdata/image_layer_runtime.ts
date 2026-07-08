import { useChannelSelectionStats } from "@spatialdata/avivatorish";
import type {
    LayerChannelSelection,
    UseLayerChannelStateResult,
} from "@spatialdata/avivatorish";
import { useLayoutEffect } from "react";

import type { StoreApi } from "zustand";
import type { ChannelsState, ViewerState } from "@/react/components/avivatorish/state";

export const DEFAULT_TONE = 0.5;

type Range = [number, number];
type RasterSlice = { width: number; height: number; data: ArrayLike<number> };

const EMPTY_RASTER: RasterSlice = { width: 0, height: 0, data: new Float32Array() };

function padTone(raw: unknown, count: number): number[] {
    return Array.from({ length: count }, (_, index) =>
        Array.isArray(raw) && typeof raw[index] === "number" ? raw[index] : DEFAULT_TONE,
    );
}

/**
 * Build padded tone arrays from the raw `brightness` / `contrast` values.
 *
 * Takes the arrays directly (not the enclosing `vivLayerProps` object) so a
 * caller in render can key off them: `vivLayerProps` is patched **in place**, so
 * its object ref is stable while the inner arrays are replaced on each tone edit.
 * Depending on the inner arrays gives the caller's `useMemo` an honest dependency
 * that recomputes exactly when tone changes. (The image panel is observer-wrapped,
 * so the React Compiler does not memoize it — the caller memoizes by hand.)
 */
export function toneFromArrays(
    brightnessRaw: unknown,
    contrastRaw: unknown,
    count: number,
): { brightness: number[]; contrast: number[] } {
    return { brightness: padTone(brightnessRaw, count), contrast: padTone(contrastRaw, count) };
}

/** Convenience wrapper that reads tone off a `vivLayerProps` object (for event handlers). */
export function toneFromVivLayerProps(
    vivLayerProps: Record<string, unknown> | undefined,
    count: number,
): { brightness: number[]; contrast: number[] } {
    return toneFromArrays(vivLayerProps?.brightness, vivLayerProps?.contrast, count);
}

function channelOptionsForSelections(
    selections: LayerChannelSelection[],
    channelNames: string[],
): string[] {
    return selections.map((selection, index) => {
        const channelIndex = selection.c ?? index;
        return channelNames[channelIndex] ?? `Channel ${channelIndex + 1}`;
    });
}

type RuntimeInput = {
    hookState: UseLayerChannelStateResult;
    loader: unknown;
    channelNames: string[];
    channelsStore: StoreApi<ChannelsState>;
    viewerStore: StoreApi<ViewerState>;
};

/**
 * One-way projection of **stats + flags** into the per-panel Avivator zustand
 * stores. The stats engine itself (load / cache / cancel, keyed by
 * channelId + selection) now lives in the library `useChannelSelectionStats`;
 * this hook only *projects* its output — histogram `domains`/`raster` + the
 * per-channel loading flag — into the stores the channel UI reads, alongside the
 * `loader`, `channelOptions` and `pixelValues` flags. Channel config
 * (colors/contrast/visibility/selections) and tone are NOT mirrored here — they
 * live in the canonical render-stack entry and the UI reads them through the
 * panel context. Never writes MobX; never subscribes zustand → persistence.
 */
export function useImageLayerRuntime({
    hookState,
    loader,
    channelNames,
    channelsStore,
    viewerStore,
}: RuntimeInput) {
    const { statsByIndex, loadingByChannelId } = useChannelSelectionStats({
        loader,
        channelIds: hookState.channelIds,
        selections: hookState.selections,
        use3d: false,
        // Show the persisted contrast limits as the histogram domain until a
        // channel's real stats land (matches the previous pre-fetch display).
        fallbackDomains: hookState.contrastLimits,
    });

    // `hookState.channelIds` / `selections` / `contrastLimits` are stable array
    // refs from `useLayerChannelState`'s store — they only change identity on an
    // actual channel edit — so they work directly as effect deps (no string-key
    // proxy needed), which keeps the dependency list honest for the linter.
    useLayoutEffect(() => {
        channelsStore.setState({
            domains: hookState.channelIds.map(
                (_, index) =>
                    statsByIndex[index]?.domain ??
                    hookState.contrastLimits[index] ??
                    ([0, 255] as Range),
            ),
            raster: hookState.channelIds.map(
                (_, index) => statsByIndex[index]?.raster ?? EMPTY_RASTER,
            ),
            loader,
        });

        const viewer = viewerStore.getState();
        viewerStore.setState({
            channelOptions: channelOptionsForSelections(hookState.selections, channelNames),
            pixelValues: hookState.channelIds.map(
                (_, index) => viewer.pixelValues[index] ?? Number.NaN,
            ),
            isChannelLoading: hookState.channelIds.map(
                (id) => loadingByChannelId.get(id) ?? false,
            ),
            isViewerLoading: false,
        });
    }, [
        channelNames,
        channelsStore,
        hookState.channelIds,
        hookState.contrastLimits,
        hookState.selections,
        loader,
        loadingByChannelId,
        statsByIndex,
        viewerStore,
    ]);
}
