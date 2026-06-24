import { Typography } from "@mui/material";
import type { LayerChannelConfig } from "@spatialdata/avivatorish";
import {
    ImageLayerContextProvider,
    useImageLayerContext,
    useLayerChannelState,
    type LayerConfig,
} from "@spatialdata/vis";
import { observer } from "mobx-react-lite";
import {
    createContext,
    useCallback,
    useContext,
    useMemo,
    useRef,
    useState,
} from "react";

import { useChart } from "@/react/context";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "@/react/components/SpatialDataMDVReact";
import {
    pickDefaultSelectionForAdd,
    useImageLayerRuntime,
} from "@/react/spatialdata/image_layer_runtime";
import { useRenderStackEntry } from "@/react/spatialdata/render_stack_control";
import { VivChannelList } from "../ColorChannelComponents";
import {
    createVivStores,
    VivProvider,
    useChannelsStoreApi,
    useViewerStoreApi,
} from "../avivatorish/state";

type ImageLayerConfig = Extract<LayerConfig, { type: "image" }>;

export type SpatialImagePanelContextValue = {
    layerId: string;
    loader: unknown;
    channelNames: string[];
    channelIds: string[];
    colors: [number, number, number][];
    contrastLimits: [number, number][];
    channelsVisible: boolean[];
    selections: Array<{ z: number; c: number; t: number }>;
    setChannels: ReturnType<typeof useLayerChannelState>["setChannels"];
    addChannel: () => void;
    addChannelWithSelection: () => void;
    removeChannel: (index: number) => void;
    patchVivLayerProps: (patch: Record<string, unknown>) => void;
    patchToneAtIndex: (index: number, key: "brightness" | "contrast", value: number) => void;
};

const SpatialImagePanelContext = createContext<SpatialImagePanelContextValue | null>(null);

export function useSpatialImagePanelContext() {
    return useContext(SpatialImagePanelContext);
}

/** @deprecated Use `useSpatialImagePanelContext` */
export const ImageLayerChannelContext = createContext<{ layerId: string } | null>(null);

/** @deprecated Use `useSpatialImagePanelContext` */
export function useImageLayerChannelContext() {
    return useContext(ImageLayerChannelContext);
}

function channelConfigKey(channels: LayerChannelConfig): string {
    return JSON.stringify(channels);
}

const ImageLayerPanelLoaded = observer(function ImageLayerPanelLoaded({
    entryId,
    elementKey,
}: {
    entryId: string;
    elementKey: string;
}) {
    const layerContext = useImageLayerContext(elementKey);
    const { layer, patchLayer } = useRenderStackEntry(entryId);

    if (!layerContext || !layer || layer.type !== "image") {
        return (
            <Typography variant="body2" color="text.secondary">
                Loading image layer…
            </Typography>
        );
    }

    return (
        <ImageLayerPanelReady
            layer={layer}
            layerContext={layerContext}
            patchLayer={patchLayer}
        />
    );
});

const ImageLayerPanelReady = observer(function ImageLayerPanelReady({
    layer,
    layerContext,
    patchLayer,
}: {
    layer: ImageLayerConfig;
    layerContext: NonNullable<ReturnType<typeof useImageLayerContext>>;
    patchLayer: ReturnType<typeof useRenderStackEntry>["patchLayer"];
}) {
    const channelsStore = useChannelsStoreApi();
    const viewerStore = useViewerStoreApi();
    const persistedKeyRef = useRef("");

    const config = layer;
    const channels = config.channels ?? {};
    const vivLayerProps =
        config.vivLayerProps && typeof config.vivLayerProps === "object"
            ? (config.vivLayerProps as Record<string, unknown>)
            : undefined;

    const onChannelsChange = useCallback(
        (next: LayerChannelConfig) => {
            const nextKey = channelConfigKey(next);
            if (nextKey === persistedKeyRef.current) return;
            persistedKeyRef.current = nextKey;
            patchLayer({ channels: next });
        },
        [patchLayer],
    );

    const hookState = useLayerChannelState({
        config: channels,
        defaults: layerContext.defaults,
        layerId: config.id,
        onChannelsChange,
    });

    useImageLayerRuntime({
        hookState,
        loader: layerContext.loader,
        channelNames: layerContext.channelNames,
        vivLayerProps,
        channelsStore,
        viewerStore,
    });

    const hookStateRef = useRef(hookState);
    hookStateRef.current = hookState;

    const patchVivLayerProps = useCallback(
        (patch: Record<string, unknown>) => {
            patchLayer({ vivLayerProps: patch });
        },
        [patchLayer],
    );

    const patchToneAtIndex = useCallback(
        (index: number, key: "brightness" | "contrast", value: number) => {
            const count = hookStateRef.current.channelIds.length;
            const { brightness: storeBrightness, contrast: storeContrast } = channelsStore.getState();
            const brightness = [...storeBrightness];
            const contrast = [...storeContrast];
            while (brightness.length < count) brightness.push(0.5);
            while (contrast.length < count) contrast.push(0.5);
            if (key === "brightness") {
                brightness[index] = value;
                patchVivLayerProps({ brightness });
            } else {
                contrast[index] = value;
                patchVivLayerProps({ contrast });
            }
        },
        [channelsStore, patchVivLayerProps],
    );

    const { addChannel, setChannels } = hookState;

    const addChannelWithSelection = useCallback(() => {
        const current = hookStateRef.current;
        const nextC = pickDefaultSelectionForAdd(
            current.selections,
            layerContext.channelNames,
        );
        const base = current.selections[0] ?? {};
        addChannel();
        setChannels({
            selections: [
                ...current.selections.map(
                    (selection: (typeof current.selections)[number]) => ({ ...selection }),
                ),
                { ...base, c: nextC },
            ],
        });
    }, [addChannel, layerContext.channelNames, setChannels]);

    const panelContext = useMemo<SpatialImagePanelContextValue>(
        () => ({
            layerId: config.id,
            loader: layerContext.loader,
            channelNames: layerContext.channelNames,
            channelIds: hookState.channelIds,
            colors: hookState.colors,
            contrastLimits: hookState.contrastLimits,
            channelsVisible: hookState.channelsVisible,
            selections: hookState.selections.map((selection: (typeof hookState.selections)[number]) => ({
                z: selection.z ?? 0,
                c: selection.c ?? 0,
                t: selection.t ?? 0,
            })),
            setChannels: hookState.setChannels,
            addChannel: hookState.addChannel,
            addChannelWithSelection,
            removeChannel: hookState.removeChannel,
            patchVivLayerProps,
            patchToneAtIndex,
        }),
        [
            addChannelWithSelection,
            config.id,
            hookState.channelIds,
            hookState.channelsVisible,
            hookState.colors,
            hookState.contrastLimits,
            hookState.addChannel,
            hookState.removeChannel,
            hookState.selections,
            hookState.setChannels,
            layerContext.channelNames,
            layerContext.loader,
            patchToneAtIndex,
            patchVivLayerProps,
        ],
    );

    return (
        <SpatialImagePanelContext.Provider value={panelContext}>
            <ImageLayerChannelContext.Provider value={{ layerId: config.id }}>
                <div className="space-y-2">
                    <Typography variant="caption" color="text.secondary">
                        Image source: SpatialData zarr
                    </Typography>
                    <VivChannelList />
                </div>
            </ImageLayerChannelContext.Provider>
        </SpatialImagePanelContext.Provider>
    );
});

const ImageLayerPanelBridge = observer(function ImageLayerPanelBridge({
    entryId,
}: {
    entryId: string;
}) {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const { entry } = useRenderStackEntry(entryId);
    const registry = chart.imageLayerRegistry;
    const elementKey = entry?.kind === "spatial" ? entry.source.elementKey : undefined;

    if (!registry || !elementKey) {
        return (
            <Typography variant="body2" color="text.secondary">
                Waiting for image data…
            </Typography>
        );
    }

    return (
        <ImageLayerContextProvider
            getImageLoadedDataByElementKey={registry.getImageLoadedDataByElementKey}
            getLayerLoadStateByElementKey={registry.getLayerLoadStateByElementKey}
        >
            <ImageLayerPanelLoaded entryId={entryId} elementKey={elementKey} />
        </ImageLayerContextProvider>
    );
});

type Props = {
    entryId: string;
};

export default function ImageLayerPanel({ entryId }: Props) {
    const [vivStores] = useState(createVivStores);
    return (
        <VivProvider vivStores={vivStores}>
            <ImageLayerPanelBridge entryId={entryId} />
        </VivProvider>
    );
}
