import { Typography } from "@mui/material";
import {
    pickDefaultSelectionForAdd,
    serializeChannelConfig,
    type LayerChannelConfig,
} from "@spatialdata/avivatorish";
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
    useRef,
    useState,
} from "react";

import { useChart } from "@/react/context";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "@/react/components/SpatialDataMDVReact";
import {
    toneFromArrays,
    toneFromVivLayerProps,
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
    /** Tone, derived from the entry's `vivLayerProps` (host-owned sibling of `channels`). */
    brightness: number[];
    contrast: number[];
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

// Stable identity for a layer with no persisted channel config yet. A fresh `{}`
// per render would change the `config` prop of `useLayerChannelState` every render,
// refiring its sync effect → setState → re-render loop (stack overflow on new layers).
const EMPTY_CHANNELS: LayerChannelConfig = Object.freeze({});

function channelConfigKey(channels: LayerChannelConfig): string {
    return serializeChannelConfig(channels);
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
    const channels = config.channels ?? EMPTY_CHANNELS;
    const vivLayerProps =
        config.vivLayerProps && typeof config.vivLayerProps === "object"
            ? (config.vivLayerProps as Record<string, unknown>)
            : undefined;
    // Read the tone arrays during render so this `observer` subscribes to the
    // in-place `vivLayerProps` patches (see render_stack_observe.touchVivLayerProps).
    // Tone is canonical here — the slider thumbs read it through the panel context.
    const brightnessRaw = vivLayerProps?.brightness;
    const contrastRaw = vivLayerProps?.contrast;

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
        channelsStore,
        viewerStore,
    });

    const channelCount = hookState.channelIds.length;
    // `vivLayerProps` is patched in place (stable object ref); `brightnessRaw` /
    // `contrastRaw` are the inner arrays, replaced on each tone edit. Computing
    // tone from them (not the object) lets the React Compiler memoize this — no
    // hand-written `useMemo`, and the in-place patch still recomputes correctly.
    const tone = toneFromArrays(brightnessRaw, contrastRaw, channelCount);

    const hookStateRef = useRef(hookState);
    hookStateRef.current = hookState;
    const vivLayerPropsRef = useRef(vivLayerProps);
    vivLayerPropsRef.current = vivLayerProps;

    const patchVivLayerProps = useCallback(
        (patch: Record<string, unknown>) => {
            patchLayer({ vivLayerProps: patch });
        },
        [patchLayer],
    );

    const patchToneAtIndex = useCallback(
        (index: number, key: "brightness" | "contrast", value: number) => {
            // Read base tone from the canonical `vivLayerProps`, not the runtime
            // store — there is a single source of truth for tone.
            const count = hookStateRef.current.channelIds.length;
            const current = toneFromVivLayerProps(vivLayerPropsRef.current, count);
            const next = key === "brightness" ? [...current.brightness] : [...current.contrast];
            next[index] = value;
            patchVivLayerProps({ [key]: next });
        },
        [patchVivLayerProps],
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

    // Plain object — the React Compiler memoizes it (and its reactive inputs)
    // without a hand-maintained dependency list.
    const panelContext: SpatialImagePanelContextValue = {
        layerId: config.id,
        loader: layerContext.loader,
        channelNames: layerContext.channelNames,
        channelIds: hookState.channelIds,
        colors: hookState.colors,
        contrastLimits: hookState.contrastLimits,
        channelsVisible: hookState.channelsVisible,
        brightness: tone.brightness,
        contrast: tone.contrast,
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
    };

    return (
        <SpatialImagePanelContext.Provider value={panelContext}>
            <div className="space-y-2">
                <Typography variant="caption" color="text.secondary">
                    Image source: SpatialData zarr
                </Typography>
                <VivChannelList />
            </div>
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
