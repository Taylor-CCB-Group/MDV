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

export function toneFromVivLayerProps(
    vivLayerProps: Record<string, unknown> | undefined,
    count: number,
): { brightness: number[]; contrast: number[] } {
    const brightnessRaw = vivLayerProps?.brightness;
    const contrastRaw = vivLayerProps?.contrast;
    return {
        brightness: Array.from({ length: count }, (_, index) => {
            if (Array.isArray(brightnessRaw) && typeof brightnessRaw[index] === "number") {
                return brightnessRaw[index];
            }
            return DEFAULT_TONE;
        }),
        contrast: Array.from({ length: count }, (_, index) => {
            if (Array.isArray(contrastRaw) && typeof contrastRaw[index] === "number") {
                return contrastRaw[index];
            }
            return DEFAULT_TONE;
        }),
    };
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

    const channelIdsKey = hookState.channelIds.join("\0");
    const selectionsKey = JSON.stringify(hookState.selections);
    const contrastLimitsKey = JSON.stringify(hookState.contrastLimits);

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
        channelIdsKey,
        channelNames,
        channelsStore,
        contrastLimitsKey,
        loader,
        loadingByChannelId,
        selectionsKey,
        statsByIndex,
        viewerStore,
    ]);
}
