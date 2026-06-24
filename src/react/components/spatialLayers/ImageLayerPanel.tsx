import { Typography } from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";
import { observer } from "mobx-react-lite";
import { createContext, useContext, useEffect, useMemo, useRef, useState } from "react";

import { VivChannelList } from "../ColorChannelComponents";
import {
    buildChannelDefaultsFromConfig,
    channelStateKey,
    mergeHydrateChannelsState,
    syncViewerChannelArraysFromStore,
    serializeChannelsFromStore,
    type ChannelConfig,
    type LoaderDefaults,
} from "@/react/spatialdata/image_layer_channel_bridge";
import { useRenderStackEntry } from "@/react/spatialdata/render_stack_control";
import {
    createVivStores,
    VivProvider,
    useChannelsStoreApi,
    useViewerStoreApi,
} from "../avivatorish/state";

type ImageLayerConfig = Extract<LayerConfig, { type: "image" }>;

export const ImageLayerChannelContext = createContext<{ layerId: string } | null>(null);

export function useImageLayerChannelContext() {
    return useContext(ImageLayerChannelContext);
}

type Props = {
    entryId: string;
    imageSource?: "spatial" | "ome_tiff";
    loaderDefaults?: LoaderDefaults;
    channelNames?: string[];
};

const ImageLayerChannelBridge = observer(function ImageLayerChannelBridge({
    entryId,
    imageSource = "spatial",
    loaderDefaults,
    channelNames = [],
}: Props) {
    const { layer, patchLayer } = useRenderStackEntry(entryId);
    const channelsStore = useChannelsStoreApi();
    const viewerStore = useViewerStoreApi();
    const hydratingRef = useRef(false);
    const persistedKeyRef = useRef("");
    const channelsRef = useRef<ChannelConfig>({});

    if (!layer || layer.type !== "image") {
        return null;
    }

    const config = layer as ImageLayerConfig;
    const channels = config.channels ?? {};
    channelsRef.current = channels;

    const defaults = useMemo(
        () =>
            buildChannelDefaultsFromConfig({
                layerId: config.id,
                channels,
                loaderDefaults,
                channelNames,
            }),
        [channels, config.id, loaderDefaults, channelNames],
    );

    useEffect(() => {
        return channelsStore.subscribe((state) => {
            if (hydratingRef.current) return;
            const nextChannels = serializeChannelsFromStore(state, channelsRef.current);
            const nextKey = channelStateKey(nextChannels);
            if (nextKey === persistedKeyRef.current) return;
            persistedKeyRef.current = nextKey;
            channelsRef.current = nextChannels;

            const viewer = viewerStore.getState();
            viewerStore.setState({
                ...syncViewerChannelArraysFromStore(state, viewer, channelNames),
                isViewerLoading: false,
            });
            patchLayer({ channels: nextChannels });
        });
    }, [channelNames, channelsStore, patchLayer, viewerStore]);

    useEffect(() => {
        const nextChannels: ChannelConfig = {
            ...channels,
            channelIds: defaults.ids,
            colors: defaults.colors,
            contrastLimits: defaults.contrastLimits,
            channelsVisible: defaults.channelsVisible,
            selections: defaults.selections.map((selection) => ({
                z: selection.z,
                c: selection.c,
                t: selection.t,
            })),
        };
        const nextKey = channelStateKey(nextChannels);
        if (nextKey === persistedKeyRef.current) return;

        const existing = channelsStore.getState();
        const expanding = defaults.ids.length > existing.ids.length;

        hydratingRef.current = true;
        persistedKeyRef.current = nextKey;
        channelsRef.current = nextChannels;

        const merged = mergeHydrateChannelsState(defaults, existing);
        const fullReplace =
            expanding ||
            defaults.ids.length !== existing.ids.length ||
            defaults.ids.some((id, index) => id !== existing.ids[index]);
        channelsStore.setState(merged);
        const viewer = viewerStore.getState();
        const viewerArrays = syncViewerChannelArraysFromStore(
            {
                ids: merged.ids,
                colors: merged.colors,
                contrastLimits: merged.contrastLimits,
                channelsVisible: merged.channelsVisible,
                selections: merged.selections,
            },
            viewer,
            channelNames,
        );
        viewerStore.setState({
            ...viewerArrays,
            ...(fullReplace
                ? {
                      isChannelLoading: merged.ids.map(() => false),
                      pixelValues: merged.ids.map(() => Number.NaN),
                  }
                : {}),
            isViewerLoading: false,
        });
        hydratingRef.current = false;
    }, [channelNames, channels, channelsStore, defaults, viewerStore]);

    return (
        <ImageLayerChannelContext.Provider value={{ layerId: config.id }}>
            <div className="space-y-2">
                <Typography variant="caption" color="text.secondary">
                    Image source: {imageSource === "spatial" ? "SpatialData zarr" : "OME-TIFF"}
                </Typography>
                <VivChannelList />
            </div>
        </ImageLayerChannelContext.Provider>
    );
});

export default function ImageLayerPanel(props: Props) {
    const [vivStores] = useState(createVivStores);
    return (
        <VivProvider vivStores={vivStores}>
            <ImageLayerChannelBridge {...props} />
        </VivProvider>
    );
}
