import { Typography } from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";
import { useEffect, useMemo, useRef, useState } from "react";

import { VivChannelList } from "../ColorChannelComponents";
import {
    createVivStores,
    type ChannelsState,
    VivProvider,
    useChannelsStoreApi,
    useViewerStoreApi,
} from "../avivatorish/state";

type ImageLayerConfig = Extract<LayerConfig, { type: "image" }>;
type ChannelConfig = NonNullable<ImageLayerConfig["channels"]>;
type SpatialSelection = NonNullable<ChannelConfig["selections"]>[number];
type VivSelection = ChannelsState["selections"][number];
type Range = [number, number];
type Rgb = [number, number, number];

type LoaderDefaults = {
    colors?: [number, number, number][];
    contrastLimits?: [number, number][];
    channelsVisible?: boolean[];
    selections?: ChannelConfig["selections"];
};

type Props = {
    config: ImageLayerConfig;
    imageSource?: "spatial" | "ome_tiff";
    loaderDefaults?: LoaderDefaults;
    channelNames?: string[];
    updateLayer: (updates: Partial<LayerConfig>) => void;
};

function fillArray<T>(values: T[] | undefined, count: number, fallback: (index: number) => T): T[] {
    return Array.from({ length: count }, (_, index) => values?.[index] ?? fallback(index));
}

function toVivSelection(selection: SpatialSelection | undefined, index: number): VivSelection {
    return {
        z: selection?.z ?? 0,
        c: selection?.c ?? index,
        t: selection?.t ?? 0,
    };
}

function fromVivSelection(selection: VivSelection): SpatialSelection {
    return {
        z: selection.z,
        c: selection.c,
        t: selection.t,
    };
}

function channelCount(
    channels: ChannelConfig,
    loaderDefaults: LoaderDefaults | undefined,
    channelNames: string[],
) {
    return Math.max(
        channels.channelIds?.length ?? 0,
        channels.colors?.length ?? 0,
        channels.contrastLimits?.length ?? 0,
        channels.channelsVisible?.length ?? 0,
        channels.selections?.length ?? 0,
        loaderDefaults?.colors?.length ?? 0,
        loaderDefaults?.contrastLimits?.length ?? 0,
        loaderDefaults?.channelsVisible?.length ?? 0,
        loaderDefaults?.selections?.length ?? 0,
        channelNames.length,
        1,
    );
}

function channelStateKey(channels: ChannelConfig) {
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

function ImageLayerChannelBridge({
    config,
    imageSource = "spatial",
    loaderDefaults,
    channelNames = [],
    updateLayer,
}: Props) {
    const channelsStore = useChannelsStoreApi();
    const viewerStore = useViewerStoreApi();
    const hydratingRef = useRef(false);
    const persistedKeyRef = useRef("");
    const channels = config.channels ?? {};
    const count = channelCount(channels, loaderDefaults, channelNames);

    const defaults = useMemo(() => {
        const colors = fillArray<Rgb>(
            channels.colors,
            count,
            (index) => loaderDefaults?.colors?.[index] ?? [255, 255, 255],
        );
        const contrastLimits = fillArray(
            channels.contrastLimits,
            count,
            (index) => loaderDefaults?.contrastLimits?.[index] ?? [0, 255] as Range,
        );
        const domains = fillArray(
            loaderDefaults?.contrastLimits,
            count,
            (index) => fallbackDomain(contrastLimits[index]),
        );
        const selections = Array.from({ length: count }, (_, index) =>
            toVivSelection(channels.selections?.[index] ?? loaderDefaults?.selections?.[index], index),
        );
        return {
            ids: fillArray(channels.channelIds, count, (index) => `${config.id}-channel-${index}`),
            colors,
            contrastLimits,
            domains,
            channelsVisible: fillArray(
                channels.channelsVisible,
                count,
                (index) => loaderDefaults?.channelsVisible?.[index] ?? true,
            ),
            selections,
            raster: fillArray(undefined, count, () => ({
                width: 0,
                height: 0,
                data: new Float32Array(),
            })),
            channelOptions: fillArray(channelNames, count, (index) => `Channel ${index + 1}`),
        };
    }, [channels, config.id, count, loaderDefaults, channelNames]);

    useEffect(() => {
        const nextChannels: ChannelConfig = {
            ...channels,
            channelIds: defaults.ids,
            colors: defaults.colors,
            contrastLimits: defaults.contrastLimits,
            channelsVisible: defaults.channelsVisible,
            selections: defaults.selections.map(fromVivSelection),
        };
        const nextKey = channelStateKey(nextChannels);
        if (nextKey === persistedKeyRef.current) return;
        hydratingRef.current = true;
        persistedKeyRef.current = nextKey;
        channelsStore.setState({
            ids: defaults.ids,
            colors: defaults.colors,
            contrastLimits: defaults.contrastLimits,
            domains: defaults.domains,
            channelsVisible: defaults.channelsVisible,
            selections: defaults.selections,
            raster: defaults.raster,
        });
        viewerStore.setState({
            channelOptions: defaults.channelOptions,
            isChannelLoading: defaults.ids.map(() => false),
            pixelValues: defaults.ids.map(() => Number.NaN),
            isViewerLoading: false,
        });
        hydratingRef.current = false;
    }, [channels, channelsStore, defaults, viewerStore]);

    useEffect(() => {
        return channelsStore.subscribe((state) => {
            if (hydratingRef.current) return;
            const nextChannels: ChannelConfig = {
                ...channels,
                channelIds: state.ids,
                colors: state.colors,
                contrastLimits: state.contrastLimits,
                channelsVisible: state.channelsVisible,
                selections: state.selections.map(fromVivSelection),
            };
            const nextKey = channelStateKey(nextChannels);
            if (nextKey === persistedKeyRef.current) return;
            persistedKeyRef.current = nextKey;
            updateLayer({ channels: nextChannels });
        });
    }, [channels, channelsStore, updateLayer]);

    return (
        <div className="space-y-2">
            <Typography variant="caption" color="text.secondary">
                Image source: {imageSource === "spatial" ? "SpatialData zarr" : "OME-TIFF"}
            </Typography>
            <VivChannelList />
        </div>
    );
}

export default function ImageLayerPanel(props: Props) {
    const [vivStores] = useState(createVivStores);
    return (
        <VivProvider vivStores={vivStores}>
            <ImageLayerChannelBridge {...props} />
        </VivProvider>
    );
}
