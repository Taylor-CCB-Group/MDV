import { getSingleSelectionStats } from "@spatialdata/avivatorish";
import type {
    LayerChannelSelection,
    UseLayerChannelStateResult,
} from "@spatialdata/avivatorish";
import { useEffect, useLayoutEffect, useRef } from "react";

import type { StoreApi } from "zustand";
import type { ChannelsState, ViewerState } from "@/react/components/avivatorish/state";

export const DEFAULT_TONE = 0.5;

type Range = [number, number];
type RasterSlice = { width: number; height: number; data: ArrayLike<number> };

function toRasterSlice(raster: { width: number; height: number; data: ArrayLike<number> } | undefined): RasterSlice {
    if (!raster) return EMPTY_RASTER;
    const data =
        raster.data instanceof Float32Array ? raster.data : Float32Array.from(raster.data);
    return { width: raster.width, height: raster.height, data };
}

type ChannelStats = {
    domains: Range;
    contrastLimits: Range;
    raster: RasterSlice;
};

const EMPTY_RASTER: RasterSlice = { width: 0, height: 0, data: new Float32Array() };

export function pickDefaultSelectionForAdd(
    selections: LayerChannelSelection[],
    channelNames: string[],
): number {
    const used = new Set(selections.map((selection) => selection.c ?? 0));
    for (let c = 0; c < channelNames.length; c++) {
        if (!used.has(c)) return c;
    }
    return 0;
}

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

/** Stats are keyed by channel row + z/c/t selection, not channelId alone. */
export function selectionStatsKey(
    channelId: string,
    selection: LayerChannelSelection | undefined,
    index: number,
): string {
    return `${channelId}:${selection?.z ?? 0}:${selection?.c ?? index}:${selection?.t ?? 0}`;
}

function buildSelectionStatsKeys(
    channelIds: string[],
    selections: LayerChannelSelection[],
): string[] {
    return channelIds.map((channelId, index) =>
        selectionStatsKey(channelId, selections[index], index),
    );
}

type RuntimeInput = {
    hookState: UseLayerChannelStateResult;
    loader: unknown;
    channelNames: string[];
    vivLayerProps?: Record<string, unknown>;
    channelsStore: StoreApi<ChannelsState>;
    viewerStore: StoreApi<ViewerState>;
};

/**
 * One-way projection from `useLayerChannelState` + stats cache into per-panel Avivator zustand stores.
 * Never writes MobX; never subscribes zustand → persistence.
 */
export function useImageLayerRuntime({
    hookState,
    loader,
    channelNames,
    vivLayerProps,
    channelsStore,
    viewerStore,
}: RuntimeInput) {
    const statsCacheRef = useRef(new Map<string, ChannelStats>());
    /** Per-row selection key last successfully loaded (not merely attempted). */
    const completedSelectionKeysRef = useRef<string[]>([]);

    const channelIdsKey = hookState.channelIds.join("\0");
    const colorsKey = JSON.stringify(hookState.colors);
    const channelsVisibleKey = JSON.stringify(hookState.channelsVisible);
    const selectionsKey = JSON.stringify(hookState.selections);
    const contrastLimitsKey = JSON.stringify(hookState.contrastLimits);
    const vivToneKey = JSON.stringify({
        brightness: vivLayerProps?.brightness,
        contrast: vivLayerProps?.contrast,
    });

    const selectionSignature = buildSelectionStatsKeys(
        hookState.channelIds,
        hookState.selections,
    ).join("\0");
    const channelCount = hookState.channelIds.length;

    useLayoutEffect(() => {
        const count = hookState.channelIds.length;
        const tone = toneFromVivLayerProps(vivLayerProps, count);
        const statsCache = statsCacheRef.current;
        const selectionKeys = buildSelectionStatsKeys(
            hookState.channelIds,
            hookState.selections,
        );

        channelsStore.setState({
            ids: [...hookState.channelIds],
            colors: hookState.colors.map((color) => [...color] as [number, number, number]),
            contrastLimits: hookState.contrastLimits.map((limits) => [...limits] as Range),
            channelsVisible: [...hookState.channelsVisible],
            selections: hookState.selections.map((selection) => ({
                z: selection.z ?? 0,
                c: selection.c ?? 0,
                t: selection.t ?? 0,
            })),
            brightness: tone.brightness,
            contrast: tone.contrast,
            domains: selectionKeys.map((key, index) => {
                const cached = statsCache.get(key);
                return cached?.domains ?? hookState.contrastLimits[index] ?? ([0, 255] as Range);
            }),
            raster: selectionKeys.map((key) => {
                const cached = statsCache.get(key);
                return cached?.raster ?? EMPTY_RASTER;
            }),
            loader,
        });

        const viewer = viewerStore.getState();
        viewerStore.setState({
            channelOptions: channelOptionsForSelections(hookState.selections, channelNames),
            pixelValues: hookState.channelIds.map(
                (_, index) => viewer.pixelValues[index] ?? Number.NaN,
            ),
            isViewerLoading: false,
        });
    }, [
        channelIdsKey,
        channelNames,
        channelsStore,
        channelsVisibleKey,
        colorsKey,
        contrastLimitsKey,
        loader,
        selectionsKey,
        viewerStore,
        vivToneKey,
    ]);

    useEffect(() => {
        if (!loader) return;

        const nextKeys = buildSelectionStatsKeys(hookState.channelIds, hookState.selections);
        if (completedSelectionKeysRef.current.length > nextKeys.length) {
            completedSelectionKeysRef.current = completedSelectionKeysRef.current.slice(
                0,
                nextKeys.length,
            );
        }

        const indicesToLoad: number[] = [];
        for (let index = 0; index < nextKeys.length; index++) {
            if (nextKeys[index] !== completedSelectionKeysRef.current[index]) {
                indicesToLoad.push(index);
            }
        }

        if (indicesToLoad.length === 0) {
            viewerStore.setState((state) => {
                const isChannelLoading = hookState.channelIds.map(() => false);
                if (
                    isChannelLoading.length === state.isChannelLoading.length &&
                    isChannelLoading.every((loading, i) => loading === state.isChannelLoading[i])
                ) {
                    return state;
                }
                return { isChannelLoading };
            });
            return;
        }

        let cancelled = false;
        const statsCache = statsCacheRef.current;
        const count = hookState.channelIds.length;

        viewerStore.setState((state) => {
            const isChannelLoading = hookState.channelIds.map(
                (_, index) =>
                    indicesToLoad.includes(index) ? true : (state.isChannelLoading[index] ?? false),
            );
            return { isChannelLoading };
        });

        const markCompleted = (index: number, key: string) => {
            const completed = [...completedSelectionKeysRef.current];
            while (completed.length <= index) completed.push("");
            completed[index] = key;
            completedSelectionKeysRef.current = completed;
        };

        const applyCachedStats = (index: number, key: string) => {
            const cached = statsCache.get(key);
            if (!cached) return false;
            channelsStore.setState((state) => {
                const domains = [...state.domains];
                const raster = [...state.raster];
                domains[index] = cached.domains;
                raster[index] = cached.raster;
                return { domains, raster };
            });
            viewerStore.setState((state) => {
                const isChannelLoading = [...state.isChannelLoading];
                while (isChannelLoading.length < count) isChannelLoading.push(false);
                isChannelLoading[index] = false;
                return { isChannelLoading };
            });
            markCompleted(index, key);
            return true;
        };

        void (async () => {
            for (const index of indicesToLoad) {
                if (cancelled) return;

                const channelId = hookState.channelIds[index];
                const selection = hookState.selections[index];
                if (!channelId) continue;

                const statsKey = nextKeys[index];
                if (applyCachedStats(index, statsKey)) continue;

                const vivSelection = {
                    z: selection?.z ?? 0,
                    c: selection?.c ?? index,
                    t: selection?.t ?? 0,
                };

                try {
                    const stats = await getSingleSelectionStats({
                        loader,
                        selection: vivSelection,
                        use3d: false,
                        includeRaster: true,
                    });
                    if (cancelled) return;

                    statsCache.set(statsKey, {
                        domains: stats.domain,
                        contrastLimits: stats.contrastLimits,
                        raster: toRasterSlice(stats.raster),
                    });

                    channelsStore.setState((state) => {
                        const domains = [...state.domains];
                        const raster = [...state.raster];
                        domains[index] = stats.domain;
                        raster[index] = toRasterSlice(stats.raster);
                        return { domains, raster };
                    });

                    viewerStore.setState((state) => {
                        const isChannelLoading = [...state.isChannelLoading];
                        while (isChannelLoading.length < count) isChannelLoading.push(false);
                        isChannelLoading[index] = false;
                        return { isChannelLoading };
                    });
                    markCompleted(index, statsKey);
                } catch (error) {
                    console.error("failed to load channel stats", error);
                    if (cancelled) return;
                    viewerStore.setState((state) => {
                        const isChannelLoading = [...state.isChannelLoading];
                        while (isChannelLoading.length < count) isChannelLoading.push(false);
                        isChannelLoading[index] = false;
                        return { isChannelLoading };
                    });
                }
            }
        })();

        return () => {
            cancelled = true;
        };
    }, [channelCount, channelsStore, loader, selectionSignature, viewerStore]);
}
