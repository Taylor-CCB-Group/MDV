import type { LayerConfig } from "@spatialdata/vis";

type ImageLayerConfig = Extract<LayerConfig, { type: "image" }>;
export type ChannelConfig = NonNullable<ImageLayerConfig["channels"]>;
type SpatialSelection = NonNullable<ChannelConfig["selections"]>[number];

export type LoaderDefaults = {
    colors?: [number, number, number][];
    contrastLimits?: [number, number][];
    channelsVisible?: boolean[];
    selections?: ChannelConfig["selections"];
};

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

const EMPTY_RASTER = { width: 0, height: 0, data: new Float32Array() };

function parseOmeroColor(color: string | undefined): Rgb | undefined {
    if (!color) return undefined;
    const hex = color.startsWith("#") ? color.slice(1) : color;
    if (hex.length !== 6) return undefined;
    const r = Number.parseInt(hex.slice(0, 2), 16);
    const g = Number.parseInt(hex.slice(2, 4), 16);
    const b = Number.parseInt(hex.slice(4, 6), 16);
    if ([r, g, b].some((value) => Number.isNaN(value))) return undefined;
    return [r, g, b];
}

export function createStableChannelId(layerId: string): string {
    const suffix =
        typeof crypto !== "undefined" && typeof crypto.randomUUID === "function"
            ? crypto.randomUUID().slice(0, 8)
            : Math.random().toString(36).slice(2, 10);
    return `${layerId}-ch-${suffix}`;
}

export function persistedChannelCount(
    channels: ChannelConfig,
    loaderDefaults: LoaderDefaults | undefined,
    channelNames: string[],
): number {
    const channelIdCount = channels.channelIds?.length ?? 0;
    if (channelIdCount > 0) return channelIdCount;

    const persistedLengths = [
        channels.colors?.length ?? 0,
        channels.contrastLimits?.length ?? 0,
        channels.channelsVisible?.length ?? 0,
        channels.selections?.length ?? 0,
    ];
    const persistedMax = Math.max(...persistedLengths, 0);
    if (persistedMax > 0) return persistedMax;

    return Math.max(
        loaderDefaults?.colors?.length ?? 0,
        loaderDefaults?.contrastLimits?.length ?? 0,
        loaderDefaults?.channelsVisible?.length ?? 0,
        loaderDefaults?.selections?.length ?? 0,
        channelNames.length,
        1,
    );
}

export function channelStateKey(channels: ChannelConfig): string {
    return JSON.stringify({
        channelIds: channels.channelIds,
        colors: channels.colors,
        contrastLimits: channels.contrastLimits,
        channelsVisible: channels.channelsVisible,
        selections: channels.selections,
    });
}

function fallbackDomain(limits: Range | undefined): Range {
    if (!limits) return [0, 255];
    return [Math.min(0, limits[0]), Math.max(255, limits[1])];
}

export function toVivSelection(selection: SpatialSelection | undefined, index: number): VivSelection {
    return {
        z: selection?.z ?? 0,
        c: selection?.c ?? index,
        t: selection?.t ?? 0,
    };
}

export function fromVivSelection(selection: VivSelection): SpatialSelection {
    return {
        z: selection.z,
        c: selection.c,
        t: selection.t,
    };
}

export type ChannelBridgeDefaults = ChannelStoreSlice & {
    domains: Range[];
    raster: { width: number; height: number; data: Float32Array }[];
    channelOptions: string[];
};

export function buildChannelDefaultsFromConfig({
    layerId,
    channels,
    loaderDefaults,
    channelNames,
}: {
    layerId: string;
    channels: ChannelConfig;
    loaderDefaults: LoaderDefaults | undefined;
    channelNames: string[];
}): ChannelBridgeDefaults {
    const count = persistedChannelCount(channels, loaderDefaults, channelNames);

    const ids = Array.from({ length: count }, (_, index) => {
        const existing = channels.channelIds?.[index];
        return existing ?? createStableChannelId(layerId);
    });

    const colors = Array.from({ length: count }, (_, index) =>
        channels.colors?.[index] ?? loaderDefaults?.colors?.[index] ?? ([255, 255, 255] as Rgb),
    );
    const contrastLimits = Array.from({ length: count }, (_, index) =>
        channels.contrastLimits?.[index] ??
            loaderDefaults?.contrastLimits?.[index] ??
            ([0, 255] as Range),
    );
    const domains = Array.from({ length: count }, (_, index) =>
        fallbackDomain(contrastLimits[index]),
    );
    const channelsVisible = Array.from({ length: count }, (_, index) =>
        channels.channelsVisible?.[index] ?? loaderDefaults?.channelsVisible?.[index] ?? true,
    );
    const selections = Array.from({ length: count }, (_, index) =>
        toVivSelection(channels.selections?.[index] ?? loaderDefaults?.selections?.[index], index),
    );

    return {
        ids,
        colors,
        contrastLimits,
        domains,
        channelsVisible,
        selections,
        raster: Array.from({ length: count }, () => ({
            width: 0,
            height: 0,
            data: new Float32Array(),
        })),
        channelOptions: channelOptionsForSelections(selections, channelNames),
    };
}

export function serializeChannelsFromStore(
    state: ChannelStoreSlice,
    persistedChannels: ChannelConfig,
): ChannelConfig {
    return {
        ...persistedChannels,
        channelIds: [...state.ids],
        colors: [...state.colors],
        contrastLimits: [...state.contrastLimits],
        channelsVisible: [...state.channelsVisible],
        selections: state.selections.map(fromVivSelection),
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

export function resizeViewerChannelArrays(
    current: ViewerChannelArrays,
    count: number,
    channelNames: string[],
): ViewerChannelArrays {
    return {
        channelOptions: Array.from(
            { length: count },
            (_, index) => current.channelOptions[index] ?? channelNames[index] ?? `Channel ${index + 1}`,
        ),
        isChannelLoading: Array.from(
            { length: count },
            (_, index) => current.isChannelLoading[index] ?? false,
        ),
        pixelValues: Array.from(
            { length: count },
            (_, index) => current.pixelValues[index] ?? Number.NaN,
        ),
    };
}

export function mergeHydrateChannelsState(
    defaults: ChannelBridgeDefaults,
    existing: {
        ids: string[];
        domains: Range[];
        raster: { width: number; height: number; data: ArrayLike<number> }[];
    },
): Pick<ChannelBridgeDefaults, "ids" | "colors" | "contrastLimits" | "channelsVisible" | "selections" | "domains" | "raster"> {
    const preserveRuntime =
        defaults.ids.length > 0 &&
        defaults.ids.length === existing.ids.length &&
        defaults.ids.every((id, index) => id === existing.ids[index]);

    if (!preserveRuntime) {
        return {
            ids: defaults.ids,
            colors: defaults.colors,
            contrastLimits: defaults.contrastLimits,
            channelsVisible: defaults.channelsVisible,
            selections: defaults.selections,
            domains: defaults.domains,
            raster: defaults.raster,
        };
    }

    const rasterById = new Map(
        existing.ids.map((id, index) => [id, existing.raster[index]] as const),
    );
    const domainsById = new Map(
        existing.ids.map((id, index) => [id, existing.domains[index]] as const),
    );

    return {
        ids: defaults.ids,
        colors: defaults.colors,
        contrastLimits: defaults.contrastLimits,
        channelsVisible: defaults.channelsVisible,
        selections: defaults.selections,
        domains: defaults.ids.map((id, index) => domainsById.get(id) ?? defaults.domains[index] ?? ([0, 255] as Range)),
        raster: defaults.ids.map((id) => {
            const preserved = rasterById.get(id);
            if (!preserved) return EMPTY_RASTER;
            return {
                width: preserved.width,
                height: preserved.height,
                data:
                    preserved.data instanceof Float32Array
                        ? preserved.data
                        : Float32Array.from(preserved.data),
            };
        }),
    };
}

export function loaderDefaultsFromImageChannels(
    channels: { label?: string; color?: string; active?: boolean; window?: { start: number; end: number } }[],
): LoaderDefaults | undefined {
    if (!channels.length) return undefined;
    return {
        colors: channels.map((channel) => parseOmeroColor(channel.color) ?? [255, 255, 255]),
        contrastLimits: channels.map((channel) => {
            const window = channel.window;
            if (window) return [window.start, window.end] as Range;
            return [0, 255] as Range;
        }),
        channelsVisible: channels.map((channel) => channel.active ?? true),
    };
}

export function channelNamesFromImageChannels(
    channels: { label?: string }[],
): string[] {
    return channels.map((channel, index) => channel.label ?? `Channel ${index + 1}`);
}
