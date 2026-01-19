import { type PropsWithChildren, createContext, useContext } from "react";
import type { loadBioformatsZarr, loadOmeTiff, loadOmeZarr } from "@vivjs-experimental/viv";
import { createStore } from "zustand";
import { useStoreWithEqualityFn } from "zustand/traditional";
import { observer } from "mobx-react-lite";
import type { EqFn, Selector, ZustandStore } from "./zustandTypes";

// what about loadOmeZarr, loadBioformatsZarr...
// ... not to mention HTJ2K?
export type OME_TIFF = Awaited<ReturnType<typeof loadOmeTiff>>;
export type OME_ZARR = Awaited<ReturnType<typeof loadOmeZarr>>;
export type BIO_ZARR = Awaited<ReturnType<typeof loadBioformatsZarr>>;
export type PixelSource = OME_TIFF | OME_ZARR | BIO_ZARR;

// --- copied straight from Avivator's code::: with notes / changes for MDV ---
import { RENDERING_MODES } from "@vivjs-experimental/viv";
import { getEntries } from "@/lib/utils";

const capitalize = (string: string) => string.charAt(0).toUpperCase() + string.slice(1);

// typing for generateToggles... not the most useful ones to have, but it's a start.
export type BooleanKeys<T> = {
    [K in keyof T]: T[K] extends boolean ? K : never;
}[keyof T];
export type TogglesReturnType<T> = {
    [K in BooleanKeys<T> as `toggle${Capitalize<string & K>}`]: () => void;
};
export type WithToggles<T> = T & TogglesReturnType<T>;
export type SetFunctionType<T> = (fn: (state: T) => T) => void;

function generateToggles<TDefaults extends {}, TState extends TDefaults>(
    defaults: TDefaults,
    // I don't know why we suddenly needed to add TState... typescript complaining otherwise
    set: SetFunctionType<TState>,
): TogglesReturnType<TDefaults> {
    const toggles: any = {};
    Object.entries(defaults).forEach(([k, v]) => {
        if (typeof v === "boolean") {
            toggles[`toggle${capitalize(k)}`] = () =>
                set((state) => ({
                    ...state,
                    [k]: !state[k as keyof TState],
                }));
        }
    });
    return toggles;
}
export type ChannelsState = {
    channelsVisible: boolean[];
    contrastLimits: [number, number][];
    colors: [number, number, number][];
    domains: [number, number][];
    // not all of these properties are always there - not all images have z/t
    selections: { z: number; c: number; t: number }[];
    loader: any; //TBD
    image: number;
    ids: string[];
    // props for VivContrastExtension
    brightness: number[];
    contrast: number[];
};

export const DEFAUlT_CHANNEL_STATE: ChannelsState = {
    channelsVisible: [] as boolean[],
    contrastLimits: [] as [number, number][],
    brightness: [] as number[], //not ending up as mobx observables... not sure why
    contrast: [] as number[],
    colors: [],
    domains: [] as [number, number][],
    selections: [] as { z: number; c: number; t: number }[],
    ids: [],
    // not for serialization... think about this.
    loader: [{ labels: [], shape: [] }],
    image: 0,
};
const DEFAUlT_CHANNEL_VALUES = {
    channelsVisible: true,
    contrastLimits: [0, 65535],
    brightness: 0.5,
    contrast: 0.5,
    colors: [255, 255, 255],
    domains: [0, 65535],
    selections: { z: 0, c: 0, t: 0 },
    ids: "",
};
const DEFAULT_IMAGE_STATE = {
    lensSelection: 0,
    colormap: "",
    renderingMode: RENDERING_MODES.MAX_INTENSITY_PROJECTION,
    resolution: 0,
    lensEnabled: false,
    zoomLock: true,
    panLock: true,
    isOverviewOn: false,
    useFixedAxis: true,
    xSlice: null as [number, number] | null,
    ySlice: null as [number, number] | null,
    zSlice: null as [number, number] | null,
    onViewportLoad: () => {},
} as const;
export type ImageState = typeof DEFAULT_IMAGE_STATE;
const DEFAULT_VIEWER_STATE = {
    isChannelLoading: [] as boolean[],
    isViewerLoading: true,
    pixelValues: [] as number[],
    isOffsetsSnackbarOn: false,
    loaderErrorSnackbar: {
        on: false,
        message: null as string | null,
    },
    isNoImageUrlSnackbarOn: false,
    isVolumeRenderingWarningOn: false,
    useLinkedView: false,
    isControllerOn: true,
    use3d: false, //not used at the moment, but should be
    useLens: false,
    useColormap: false,
    globalSelection: { z: 0, t: 0 },
    channelOptions: [] as any[],
    /** type is WIP */
    metadata: null as Metadata | null,
    viewState: null as any,
    source: undefined as { urlOrFile: string; description: string } | undefined,
    pyramidResolution: 0,
};
export type ViewerState = typeof DEFAULT_VIEWER_STATE;

// --- following how VivViewerMDV _parseChannels() works: (not used) ---
// export type MdvVivChannelConfig = {
//   name: string,
//   selections?: { c: number, z: number, t: number }[],
//   color?: `#${string}` | [r: number, g: number, b: number],
//   visible?: boolean,
//   contrastLimits?: [min: number, max: number],
//   domains?: [min: number, max: number],
// }

export type VivConfig = {
    // settling on using the same structure as Avivator for now...
    // todo more specification of which bits we want to exclude / include
    // channels: MdvVivChannelConfig[],
    // vs...
    // image_properties: ChannelsState,
    // vs...
    viewerStore?: Partial<ViewerState>;
    channelsStore?: Partial<ChannelsState>;
    imageSettingsStore?: Partial<ImageState>;
    //>>>?
    // url: string;
};

export type ROI = {
    min_x: number;
    min_y: number;
    max_x: number;
    max_y: number;
};
export type NewChannelValues = Partial<typeof DEFAUlT_CHANNEL_VALUES>;
/**
 * In Avivator, there is a zustand store for the channels, image settings, and viewers.
 * While there can be plural 'viewers', it is always of a single image, with shared channels and image settings.
 *
 * So in order to have appropriate state for multiple MDV charts, we have a context that contains the stores for each.
 * Within this context, we should be able to use much the same API as Avivator.
 */
export type VivContextType = {
    channelsStore: ZustandStore<
        WithToggles<ChannelsState> & {
            toggleIsOn: (index: number) => void;
            setPropertiesForChannel: (
                channel: number,
                newProperties: NewChannelValues,
            ) => void;
            removeChannel: (channel: number) => void;
            addChannel: (
                newProperties: NewChannelValues,
            ) => void;
        }
    >;
    imageSettingsStore: ZustandStore<WithToggles<ImageState>>;
    viewerStore: ZustandStore<
        WithToggles<ViewerState> & {
            setIsChannelLoading: (index: number, val: boolean) => void;
            addIsChannelLoading: (val: boolean) => void;
            removeIsChannelLoading: (index: number) => void;
        }
    >; //also has some extra methods... setIsChannelLoading, addIsChannelLoading, removeIsChannelLoading
};
// ChannelsState has a bunch of arrays, we need something that describes entries in those arrays...

export function createVivStores() {
    const channelsStore = createStore((set) => ({
        ...DEFAUlT_CHANNEL_STATE,
        // ...chart.config.viv.channelsStore,
        ...generateToggles(DEFAUlT_CHANNEL_VALUES, set),
        toggleIsOn: (index: number) =>
            set((state: ChannelsState) => {
                const channelsVisible = [...state.channelsVisible];
                channelsVisible[index] = !channelsVisible[index];
                return { ...state, channelsVisible };
            }),
        setPropertiesForChannel: (channel: number, newProperties: NewChannelValues) =>
            set((state: ChannelsState) => {
                const entries = getEntries(newProperties);
                const newState: Partial<ChannelsState> = {};
                //@ts-expect-error there is a strong case for changing the structure of the state...
                entries.forEach(([property, value]) => {
                    //@ts-expect-error there is a strong case for changing the structure of the state...
                    newState[property] = [...state[property]];
                    //@ts-expect-error there is a strong case for changing the structure of the state...
                    newState[property][channel] = value;
                });
                return { ...state, ...newState };
            }),
        removeChannel: (channel: number) =>
            set((state: ChannelsState) => {
                const newState = {};
                const channelKeys = Object.keys(DEFAUlT_CHANNEL_VALUES);
                Object.keys(state).forEach((key) => {
                    if (channelKeys.includes(key)) {
                        //@ts-expect-error viv removeChannel
                        newState[key] = state[key].filter(
                            //@ts-expect-error viv removeChannel
                            (_, j) => j !== channel,
                        );
                    }
                });
                return { ...state, ...newState };
            }),
        addChannel: (newProperties: NewChannelValues) =>
            set((state: ChannelsState) => {
                const entries = getEntries(newProperties);
                const newState = { ...state };
                entries.forEach(([property, value]) => {
                    //@ts-expect-error there is a strong case for changing the structure of the state...
                    newState[property] = [...state[property], value];
                });
                Object.entries(DEFAUlT_CHANNEL_VALUES).forEach(([k, v]) => {
                    //@ts-expect-error viv addChannel
                    if (newState[k].length < newState[entries[0][0]].length) {
                        //@ts-expect-error viv addChannel
                        newState[k] = [...state[k], v];
                    }
                });
                return newState;
            }),
    })); // satisfies VivContextType['channelsStore'];
    const imageSettingsStore = createStore((set) => ({
        ...DEFAULT_IMAGE_STATE,
        // ...chart.config.viv.imageSettingsStore,
        ...generateToggles(DEFAULT_IMAGE_STATE, set),
    })) satisfies VivContextType["imageSettingsStore"];
    const viewerStore = createStore((set) => ({
        ...DEFAULT_VIEWER_STATE,
        // ...chart.config.viv.viewerStore,
        ...generateToggles(DEFAULT_VIEWER_STATE, set),
        setIsChannelLoading: (index, val) =>
            set((state) => {
                const newIsChannelLoading = [...state.isChannelLoading];
                newIsChannelLoading[index] = val;
                return { ...state, isChannelLoading: newIsChannelLoading };
            }),
        addIsChannelLoading: (val) =>
            set((state) => {
                const newIsChannelLoading = [...state.isChannelLoading, val];
                return { ...state, isChannelLoading: newIsChannelLoading };
            }),
        removeIsChannelLoading: (index) =>
            set((state) => {
                const newIsChannelLoading = [...state.isChannelLoading];
                newIsChannelLoading.splice(index, 1);
                return { ...state, isChannelLoading: newIsChannelLoading };
            }),
    })) satisfies VivContextType["viewerStore"];
    // using `as` instead of `satisfies` because of a bug in the type-checking...
    // not sure why 'channelsStore' is misbehaving when other types are ok.
    return { viewerStore, channelsStore, imageSettingsStore } as VivContextType;
}

const VivContext = createContext<VivContextType>(null as any);
/**
 * This implictly assumes that we are in a context where there is a chart, and that we
 * can access it with `useChart()`. That design decision may be revisited (and it means
 * this detail is deviating somewhat from Avivator).
 *
 * Separate providers - referring to the same `chart.vivStores` - are used both by the chart
 * itself and by the color-change GUI dialog.
 */
export const VivProvider = observer(
    ({
        children,
        vivStores,
    }: PropsWithChildren & { vivStores: VivContextType }) => {
        return (
            <VivContext.Provider value={vivStores}>
                {children}
            </VivContext.Provider>
        );
    },
);

export type StoreName = keyof VivContextType;
export type ImageSettingsStore = VivContextType["imageSettingsStore"];
export type ViewerStore = VivContextType["viewerStore"];
export type ChannelsStore = VivContextType["channelsStore"];

// I don't think it's possible to make things like useChannelsStore.setState() work...
// so we have these extra functions...
function useStoreApi<S extends StoreName>(storeName: S): VivContextType[S] {
    const store = useContext(VivContext);
    if (!store) throw "Viv state hooks must be used within a VivProvider";
    return store[storeName];
}
export const useChannelsStoreApi = () => useStoreApi("channelsStore");
export const useImageSettingsStoreApi = () => useStoreApi("imageSettingsStore");
export const useViewerStoreApi = () => useStoreApi("viewerStore");

/** should be more-or-less equivalent to equivalent avivator hook -
 * but there can be multiple viv viewers, so we have context for that.
 */
export function useChannelsStore<U>(
    selector: Selector<ChannelsStore, U>,
    equalityFn?: EqFn<U>,
) {
    const store = useChannelsStoreApi();
    return useStoreWithEqualityFn(store, selector, equalityFn);
}
export function useImageSettingsStore<U>(
    selector: Selector<ImageSettingsStore, U>,
    equalityFn?: EqFn<U>,
) {
    const store = useImageSettingsStoreApi();
    return useStoreWithEqualityFn(store, selector, equalityFn);
}

export function useViewerStore<U>(
    selector: Selector<ViewerStore, U>,
    equalityFn?: EqFn<U>,
) {
    const store = useViewerStoreApi();
    return useStoreWithEqualityFn(store, selector, equalityFn);
}

export const useLoader = () => {
    const [fullLoader, image] = useChannelsStore((store) => [
        store.loader,
        store.image,
    ]);
    return Array.isArray(fullLoader[0]) ? fullLoader[image] : fullLoader;
};
//! todo review the typing here...
type OME_METADATA = OME_ZARR['metadata'] & {
    // in practice, we seem to get something that looks like this...
    // at least, that was true for the first sample I looked at...
    // and at least that allows us to remove some ts-expect-error
    Pixels: {
        Channels: Array<{Name: string, SamplesPerPixel: number, Color?: any}>,
        //! I don't think we actually do see these on OME-ZARR
        PhysicalSizeX?: number;
        PhysicalSizeXUnit?: string;
    }
}
export type Metadata = OME_TIFF['metadata'] | OME_METADATA | BIO_ZARR['metadata'];
//export type Metadata = TiffPreviewProps["metadata"];
export const useMetadata = (): Metadata | undefined | null => {
    try {
        const image = useChannelsStore((store) => store.image);
        const metadata = useViewerStore((store) => store.metadata);
        return Array.isArray(metadata) ? metadata[image] : metadata;
    } catch (e) {
        // we now sometimes call this hook outside of a Viv context, so this is expected.
        // console.error("no metadata", e);
    }
};

/** Add default values to a *config* (not an instance of actual store) that may have been serialized before
 * certain properties were added to the model. For example, brightness and contrast did not exist,
 * there should be arrays of correct length corresponding to the number of channels.
 */
export const applyDefaultChannelState = (config: Partial<VivConfig>) => {
    // return config as VivConfig;
    const newConfig = config as VivConfig;
    if (!newConfig.channelsStore) newConfig.channelsStore = {};
    
    // Determine the expected number of channels from channelsVisible array
    const expectedLength = newConfig.channelsStore.channelsVisible?.length || 0;
    
    // Ensure ids array exists and has the correct length
    // Ignore any serialised value in config (not expected to be there other than in test config used during dev)
    if (expectedLength > 0) {
        newConfig.channelsStore.ids = Array.from({ length: expectedLength }, 
            () => String(Math.random()));
    } else {
        newConfig.channelsStore.ids = [];
    }
    
    // Handle brightness and contrast arrays (legacy support)
    for (const [k, v] of Object.entries(DEFAUlT_CHANNEL_VALUES)) {
        if (k === "ids") continue; // Already handled above
        // it would be nice if this was more generic, but for now we are explicitly dealing with brightness and contrast...
        if (k !== "brightness" && k !== "contrast") continue;
        if (newConfig.channelsStore[k] === undefined) {
            console.log(`adding default value ${v} for`, k);
            newConfig.channelsStore[k] = new Array(expectedLength).fill(v);
        }
    }
    return newConfig;
};
