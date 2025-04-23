import { useCallback, useEffect } from "react";
// import { useDropzone as useReactDropzone } from 'react-dropzone';
import { shallow } from "zustand/shallow";
// eslint-disable-next-line camelcase
import { unstable_batchedUpdates } from "react-dom";

import {
    useImageSettingsStore,
    useLoader,
    useMetadata,
    useViewerStore,
    useViewerStoreApi,
    useChannelsStoreApi,
    useImageSettingsStoreApi,
    // type PixelSource,
} from "./state";
import {
    createLoader,
    buildDefaultSelection,
    guessRgb,
    getMultiSelectionStats,
    getBoundingCube,
    isInterleaved,
} from "./utils";
import { COLOR_PALLETE, FILL_PIXEL_VALUE } from "./constants";
import { useVivConfig } from "@/react/context";

const useSavedVivConfig = () => {
    const c = useVivConfig();
    const { viewerStore, channelsStore, imageSettingsStore } = c;
    const viewerStoreApi = useViewerStoreApi();
    const channelsStoreApi = useChannelsStoreApi();
    const imageSettingsStoreApi = useImageSettingsStoreApi();
    const applyConfig = useCallback(() => {
        viewerStoreApi.setState((state) => {
            const newState = { ...state, ...viewerStore };
            return newState;
        });
        channelsStoreApi.setState((state) => {
            const newState = { ...state, ...channelsStore };
            return newState;
        });
        imageSettingsStoreApi.setState((state) => {
            const newState = { ...state, ...imageSettingsStore };
            return newState;
        });
    }, [
        viewerStore,
        channelsStore,
        imageSettingsStore,
        channelsStoreApi.setState,
        viewerStoreApi.setState,
        imageSettingsStoreApi.setState,
    ]);
    return applyConfig;
};

export const useImage = (
    source?: { description: string; urlOrFile: string },
    history?: any,
) => {
    const applyConfig = useSavedVivConfig();
    const [use3d, toggleUse3d, toggleIsOffsetsSnackbarOn] = useViewerStore(
        (store) => [
            store.use3d,
            store.toggleUse3d,
            store.toggleIsOffsetsSnackbarOn,
        ],
        shallow,
    );
    const [lensEnabled, toggleLensEnabled] = useImageSettingsStore(
        (store) => [store.lensEnabled, store.toggleLensEnabled],
        shallow,
    );
    const loader = useLoader();
    const metadata = useMetadata();
    const viewerStore = useViewerStoreApi();
    const channelsStore = useChannelsStoreApi();
    const imageSettingsStore = useImageSettingsStoreApi();
    // biome-ignore lint/correctness/useExhaustiveDependencies: disabled in viv as well, and would cause a bunch of re-running...
    useEffect(() => {
        if (!source) return;
        async function changeLoader() {
            // Placeholder
            viewerStore.setState({ isChannelLoading: [true] });
            viewerStore.setState({ isViewerLoading: true });
            if (use3d) toggleUse3d();
            if (!source) throw "this should never happen - this is a type-guard";
            const { urlOrFile } = source;
            const newLoader = await createLoader(
                urlOrFile,
                toggleIsOffsetsSnackbarOn,
                (message) =>
                    viewerStore.setState({
                        loaderErrorSnackbar: { on: true, message },
                    }),
            );
            //@ts-ignore flagging so we can get back to this
            let nextMeta: any;
            //@ts-ignore flagging so we can get back to this
            let nextLoader: any;//PixelSource | PixelSource[] | null = null;
            if (Array.isArray(newLoader)) {
                if (newLoader.length > 1) {
                    nextMeta = newLoader.map((l) => l.metadata);
                    nextLoader = newLoader.map((l) => l.data);
                } else {
                    nextMeta = newLoader[0].metadata;
                    nextLoader = newLoader[0].data;
                }
            } else {
                nextMeta = newLoader.metadata;
                nextLoader = newLoader.data;
            }
            if (nextLoader) {
                console.info(
                    "Metadata (in JSON-like form) for current file being viewed: ",
                    nextMeta,
                );
                unstable_batchedUpdates(() => {
                    channelsStore.setState({ loader: nextLoader });
                    viewerStore.setState({
                        metadata: nextMeta,
                    });
                });
                if (use3d) toggleUse3d();
                // eslint-disable-next-line no-unused-expressions
                history?.push(
                    typeof urlOrFile === "string"
                        ? `?image_url=${urlOrFile}`
                        : "",
                );
            }
        }
        if (source) changeLoader().then(applyConfig);
    }, [source, history]); // eslint-disable-line react-hooks/exhaustive-deps
    // biome-ignore lint/correctness/useExhaustiveDependencies: suppressed in viv as well, and would cause a bunch of re-running...
    useEffect(() => {
        if (!metadata) return;
        const changeSettings = async () => {
            // Placeholder
            viewerStore.setState({ isChannelLoading: [true] });
            viewerStore.setState({ isViewerLoading: true });
            if (use3d) toggleUse3d();
            const newSelections = buildDefaultSelection(loader[0]);
            const { Channels } = metadata.Pixels;
            const channelOptions = Channels.map(
                (c, i) => c.Name ?? `Channel ${i}`,
            );
            // Default RGB.
            type Limits = [number, number][];
            type Colors = [number, number, number][];
            let newContrastLimits: Limits = [];
            let newDomains: Limits = [];
            let newColors: Colors = [];
            const isRgb = guessRgb(metadata);
            if (isRgb) {
                if (isInterleaved(loader[0].shape)) {
                    // These don't matter because the data is interleaved.
                    newContrastLimits = [[0, 255]];
                    newDomains = [[0, 255]];
                    newColors = [[255, 0, 0]];
                } else {
                    newContrastLimits = [
                        [0, 255],
                        [0, 255],
                        [0, 255],
                    ];
                    newDomains = [
                        [0, 255],
                        [0, 255],
                        [0, 255],
                    ];
                    newColors = [
                        [255, 0, 0],
                        [0, 255, 0],
                        [0, 0, 255],
                    ];
                }
                if (lensEnabled) {
                    toggleLensEnabled();
                }
                viewerStore.setState({ useColormap: false, useLens: false });
            } else {
                const stats = await getMultiSelectionStats({
                    loader,
                    selections: newSelections,
                    use3d: false,
                });
                newDomains = stats.domains;
                newContrastLimits = stats.contrastLimits;
                // If there is only one channel, use white.
                newColors =
                    newDomains.length === 1
                        ? [[255, 255, 255]] as Colors
                        : newDomains.map(
                              (_, i) =>
                                  (Channels[i]?.Color?.slice(0, -1) ??
                                  COLOR_PALLETE[i]) as Colors[number],
                          );
                viewerStore.setState({
                    useLens: channelOptions.length !== 1,
                    useColormap: true,
                });
            }
            channelsStore.setState({
                ids: newDomains.map(() => String(Math.random())),
                selections: newSelections,
                domains: newDomains,
                contrastLimits: newContrastLimits,
                colors: newColors,
                channelsVisible: newColors.map(() => true),
            });
            viewerStore.setState({
                isChannelLoading: newSelections.map((i) => !i),
                isViewerLoading: false,
                pixelValues: new Array(newSelections.length).fill(
                    FILL_PIXEL_VALUE,
                ),
                // Set the global selections (needed for the UI). All selections have the same global selection.
                globalSelection: newSelections[0],
                channelOptions,
            });
            const [xSlice, ySlice, zSlice] = getBoundingCube(loader);
            imageSettingsStore.setState({
                xSlice,
                ySlice,
                zSlice,
            });
        };
        if (metadata) changeSettings().then(applyConfig);
    }, [loader, metadata]); // eslint-disable-line react-hooks/exhaustive-deps
};

// export const useDropzone = () => {
//     const viewerStore = useViewerStoreApi();
//     const handleSubmitFile = files => {
//         let newSource;
//         if (files.length === 1) {
//             newSource = {
//                 urlOrFile: files[0],
//                 // Use the trailing part of the URL (file name, presumably) as the description.
//                 description: files[0].name
//             };
//         } else {
//             newSource = {
//                 urlOrFile: files,
//                 description: 'data.zarr'
//             };
//         }
//         viewerStore.setState({ source: newSource });
//     };
//     return useReactDropzone({
//         onDrop: handleSubmitFile
//     });
// };
