import type { LayerChannelSelection } from "@spatialdata/avivatorish";

type Range = [number, number];
type Rgb = [number, number, number];

export type VivSelection = { z: number; c: number; t: number };

export type ChannelStoreSlice = {
    ids: string[];
    colors: Rgb[];
    contrastLimits: Range[];
    channelsVisible: boolean[];
    selections: VivSelection[];
};

export type ViewerChannelArrays = {
    channelOptions: string[];
    isChannelLoading: boolean[];
    pixelValues: number[];
};

export type ChannelRuntimeStats = {
    domains: Range;
    raster: { width: number; height: number; data: Float32Array };
};

export type ChannelRuntimeCache = Map<string, ChannelRuntimeStats>;

const EMPTY_RASTER = { width: 0, height: 0, data: new Float32Array() };
const DEFAULT_TONE = 0.5;

export function channelStateKey(channels: {
    channelIds?: string[];
    colors?: Rgb[];
    contrastLimits?: Range[];
    channelsVisible?: boolean[];
    selections?: LayerChannelSelection[];
}): string {
    return JSON.stringify({
        channelIds: channels.channelIds,
        colors: channels.colors,
        contrastLimits: channels.contrastLimits,
        channelsVisible: channels.channelsVisible,
        selections: channels.selections,
    });
}

export function createStableChannelId(layerId: string): string {
    const suffix =
        typeof crypto !== "undefined" && typeof crypto.randomUUID === "function"
            ? crypto.randomUUID().slice(0, 8)
            : Math.random().toString(36).slice(2, 10);
    return `${layerId}-ch-${suffix}`;
}

export function toVivSelection(selection: LayerChannelSelection | undefined, index: number): VivSelection {
    return {
        z: selection?.z ?? 0,
        c: selection?.c ?? index,
        t: selection?.t ?? 0,
    };
}

export function channelOptionsForSelections(
    selections: VivSelection[],
    channelNames: string[],
): string[] {
    return selections.map((selection, index) => {
        const channelIndex = selection.c ?? index;
        return channelNames[channelIndex] ?? `Channel ${channelIndex + 1}`;
    });
}

export function syncViewerChannelArraysFromStore(
    state: ChannelStoreSlice,
    current: ViewerChannelArrays,
    channelNames: string[],
): ViewerChannelArrays {
    const count = state.ids.length;
    const channelOptions = channelOptionsForSelections(state.selections, channelNames);
    const shrinking = count < current.channelOptions.length;
    return {
        channelOptions,
        isChannelLoading: Array.from({ length: count }, (_, index) =>
            shrinking ? false : (current.isChannelLoading[index] ?? false),
        ),
        pixelValues: Array.from({ length: count }, (_, index) =>
            shrinking ? Number.NaN : (current.pixelValues[index] ?? Number.NaN),
        ),
    };
}

export function fallbackDomain(limits: Range | undefined): Range {
    if (!limits) return [0, 255];
    return [Math.min(0, limits[0]), Math.max(255, limits[1])];
}

export function projectRuntimeCacheToArrays(
    channelIds: string[],
    cache: ChannelRuntimeCache,
    contrastLimits: Range[],
): { domains: Range[]; raster: ChannelRuntimeStats["raster"][] } {
    return {
        domains: channelIds.map((id, index) =>
            cache.get(id)?.domains ?? fallbackDomain(contrastLimits[index]),
        ),
        raster: channelIds.map((id) => cache.get(id)?.raster ?? EMPTY_RASTER),
    };
}

export function pruneRuntimeCache(cache: ChannelRuntimeCache, channelIds: string[]) {
    const keep = new Set(channelIds);
    for (const id of cache.keys()) {
        if (!keep.has(id)) {
            cache.delete(id);
        }
    }
}

export function readVivToneArrays(
    vivLayerProps: Record<string, unknown> | undefined,
    channelCount: number,
): { brightness: number[]; contrast: number[] } {
    const brightnessRaw = vivLayerProps?.brightness;
    const contrastRaw = vivLayerProps?.contrast;
    const brightness = Array.from({ length: channelCount }, (_, index) => {
        if (Array.isArray(brightnessRaw) && typeof brightnessRaw[index] === "number") {
            return brightnessRaw[index];
        }
        return DEFAULT_TONE;
    });
    const contrast = Array.from({ length: channelCount }, (_, index) => {
        if (Array.isArray(contrastRaw) && typeof contrastRaw[index] === "number") {
            return contrastRaw[index];
        }
        return DEFAULT_TONE;
    });
    return { brightness, contrast };
}

export function patchVivToneArrays(
    vivLayerProps: Record<string, unknown> | undefined,
    channelCount: number,
    index: number,
    patch: { brightness?: number; contrast?: number },
): Record<string, unknown> {
    const { brightness, contrast } = readVivToneArrays(vivLayerProps, channelCount);
    if (patch.brightness !== undefined) {
        brightness[index] = patch.brightness;
    }
    if (patch.contrast !== undefined) {
        contrast[index] = patch.contrast;
    }
    return {
        ...(vivLayerProps ?? {}),
        brightness,
        contrast,
    };
}

export function selectionStatsToRuntime(stats: {
    domain: Range;
    contrastLimits: Range;
    raster?: { width: number; height: number; data: ArrayLike<number> };
}): ChannelRuntimeStats {
    const raster = stats.raster;
    return {
        domains: stats.domain,
        raster: raster
            ? {
                  width: raster.width,
                  height: raster.height,
                  data:
                      raster.data instanceof Float32Array
                          ? raster.data
                          : Float32Array.from(raster.data),
              }
            : EMPTY_RASTER,
    };
}
