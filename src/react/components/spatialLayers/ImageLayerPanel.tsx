import { Typography } from "@mui/material";
import type { LayerChannelConfig, LayerChannelDefaults, LayerChannelSelection } from "@spatialdata/avivatorish";
import { mergeLayerChannelState, useLayerChannelState } from "@spatialdata/vis";
import type { ImageLoaderData, LayerConfig } from "@spatialdata/vis";
import { observer } from "mobx-react-lite";
import { createContext, useCallback, useContext, useEffect, useMemo, useRef, useState, type RefObject } from "react";

import { VivChannelList } from "../ColorChannelComponents";
import { useChart } from "../../context";
import {
    channelStateKey,
    patchVivToneArrays,
    pruneRuntimeCache,
    projectRuntimeCacheToArrays,
    readVivToneArrays,
    syncViewerChannelArraysFromStore,
    toVivSelection,
    type ChannelRuntimeCache,
} from "@/react/spatialdata/image_layer_runtime";
import {
    createVivStores,
    VivProvider,
    useChannelsStoreApi,
    useViewerStoreApi,
} from "../avivatorish/state";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "../SpatialDataMDVReact";

type ImageLayerConfig = Extract<LayerConfig, { type: "image" }>;

export type SpatialImagePanelContextValue = {
    isSpatialPanel: true;
    layerId: string;
    elementKey: string;
    channelNames: string[];
    setChannels: (patch: Partial<LayerChannelConfig>) => void;
    addChannel: () => void;
    removeChannel: (index: number) => void;
    updateTone: (index: number, patch: { brightness?: number; contrast?: number }) => void;
    runtimeCacheRef: RefObject<ChannelRuntimeCache>;
    bumpRuntimeProjection: () => void;
};

const SpatialImagePanelContext = createContext<SpatialImagePanelContextValue | null>(null);

export function useImageLayerChannelContext() {
    return useContext(SpatialImagePanelContext);
}

type Props = {
    config: ImageLayerConfig;
    elementKey: string;
    updateLayer: (updates: Partial<LayerConfig>) => void;
};

function useStableImageLoaderData(imageData: ImageLoaderData | undefined) {
    const stableRef = useRef<ImageLoaderData | undefined>(undefined);
    if (!imageData) {
        stableRef.current = undefined;
        return undefined;
    }
    if (stableRef.current?.loader !== imageData.loader) {
        stableRef.current = imageData;
    }
    return stableRef.current;
}

const ImageLayerPanelContent = function ImageLayerPanelContent({
    config,
    elementKey,
    imageData,
    updateLayer,
}: Props & { imageData: ImageLoaderData }) {
    const channelsStore = useChannelsStoreApi();
    const viewerStore = useViewerStoreApi();
    const runtimeCacheRef = useRef<ChannelRuntimeCache>(new Map());
    const [runtimeRevision, setRuntimeRevision] = useState(0);
    const bumpRuntimeProjection = useCallback(() => {
        setRuntimeRevision((revision) => revision + 1);
    }, []);

    const loaderDefaults = useMemo((): LayerChannelDefaults => ({
        colors: imageData.colors,
        contrastLimits: imageData.contrastLimits,
        channelsVisible: imageData.channelsVisible,
        selections: imageData.selections,
        selectionAxisSizes: imageData.selectionAxisSizes,
    }), [
        imageData.colors,
        imageData.contrastLimits,
        imageData.channelsVisible,
        imageData.selections,
        imageData.selectionAxisSizes,
    ]);

    const channelState = useLayerChannelState({
        config: config.channels ?? {},
        defaults: loaderDefaults,
        layerId: config.id,
        onChannelsChange: (next: LayerChannelConfig) => updateLayer({ channels: next }),
    });

    const mergedChannels = useMemo(
        () => mergeLayerChannelState(config.channels ?? {}, loaderDefaults, config.id),
        [config.channels, config.id, loaderDefaults],
    );

    const channelNames = imageData.channelNames ?? [];

    const updateTone = useCallback(
        (index: number, patch: { brightness?: number; contrast?: number }) => {
            const nextVivLayerProps = patchVivToneArrays(
                config.vivLayerProps,
                channelState.channelCount,
                index,
                patch,
            );
            updateLayer({ vivLayerProps: nextVivLayerProps });
            const tone = readVivToneArrays(nextVivLayerProps, channelState.channelCount);
            channelsStore.setState({
                brightness: tone.brightness,
                contrast: tone.contrast,
            });
        },
        [channelState.channelCount, channelsStore, config.vivLayerProps, updateLayer],
    );

    const spatialContext = useMemo((): SpatialImagePanelContextValue => ({
        isSpatialPanel: true,
        layerId: config.id,
        elementKey,
        channelNames,
        setChannels: channelState.setChannels,
        addChannel: channelState.addChannel,
        removeChannel: channelState.removeChannel,
        updateTone,
        runtimeCacheRef,
        bumpRuntimeProjection,
    }), [
        bumpRuntimeProjection,
        channelNames,
        channelState.addChannel,
        channelState.removeChannel,
        channelState.setChannels,
        config.id,
        elementKey,
        updateTone,
    ]);

    const channelSyncKey = useMemo(() => {
        const tone = readVivToneArrays(config.vivLayerProps, mergedChannels.channelCount);
        return JSON.stringify({
            channels: channelStateKey({
                channelIds: mergedChannels.channelIds,
                colors: mergedChannels.colors,
                contrastLimits: mergedChannels.contrastLimits,
                channelsVisible: mergedChannels.channelsVisible,
                selections: mergedChannels.selections,
            }),
            tone,
            runtimeRevision,
        });
    }, [mergedChannels, config.vivLayerProps, runtimeRevision]);

    useEffect(() => {
        pruneRuntimeCache(runtimeCacheRef.current, mergedChannels.channelIds);
        const { domains, raster } = projectRuntimeCacheToArrays(
            mergedChannels.channelIds,
            runtimeCacheRef.current,
            mergedChannels.contrastLimits,
        );
        const tone = readVivToneArrays(config.vivLayerProps, mergedChannels.channelCount);
        const selections = mergedChannels.selections.map((selection: LayerChannelSelection, index: number) =>
            toVivSelection(selection, index),
        );

        channelsStore.setState({
            ids: [...mergedChannels.channelIds],
            colors: mergedChannels.colors.map((color: [number, number, number]) => [...color] as [number, number, number]),
            contrastLimits: mergedChannels.contrastLimits.map((limits: [number, number]) => [...limits] as [number, number]),
            channelsVisible: [...mergedChannels.channelsVisible],
            selections,
            domains,
            raster,
            brightness: tone.brightness,
            contrast: tone.contrast,
            loader: imageData.loader,
        });

        const viewer = viewerStore.getState();
        const existing = channelsStore.getState();
        const fullReplace =
            mergedChannels.channelIds.length !== existing.ids.length ||
            mergedChannels.channelIds.some((id: string, index: number) => id !== existing.ids[index]);
        const viewerArrays = syncViewerChannelArraysFromStore(
            {
                ids: mergedChannels.channelIds,
                colors: mergedChannels.colors,
                contrastLimits: mergedChannels.contrastLimits,
                channelsVisible: mergedChannels.channelsVisible,
                selections,
            },
            viewer,
            channelNames,
        );
        viewerStore.setState({
            ...viewerArrays,
            ...(fullReplace
                ? {
                      isChannelLoading: mergedChannels.channelIds.map(() => false),
                      pixelValues: mergedChannels.channelIds.map(() => Number.NaN),
                  }
                : {}),
            isViewerLoading: false,
        });
    }, [channelNames, channelSyncKey, channelsStore, imageData.loader, mergedChannels, viewerStore]);

    return (
        <SpatialImagePanelContext.Provider value={spatialContext}>
            <div className="space-y-2">
                <Typography variant="caption" color="text.secondary">
                    Image source: SpatialData zarr
                </Typography>
                <VivChannelList />
            </div>
        </SpatialImagePanelContext.Provider>
    );
};

const ImageLayerPanelBody = observer(function ImageLayerPanelBody(props: Props) {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const rawImageData = chart.imageLayerRegistry?.getImageLoadedDataByElementKey(props.elementKey);
    const imageData = useStableImageLoaderData(rawImageData);

    if (!chart.imageLayerRegistry || !imageData) {
        return (
            <Typography variant="body2" color="text.secondary">
                Loading image channel data…
            </Typography>
        );
    }

    return <ImageLayerPanelContent {...props} imageData={imageData} />;
});

export default function ImageLayerPanel(props: Props) {
    const [vivStores] = useState(createVivStores);
    return (
        <VivProvider vivStores={vivStores}>
            <ImageLayerPanelBody {...props} />
        </VivProvider>
    );
}
