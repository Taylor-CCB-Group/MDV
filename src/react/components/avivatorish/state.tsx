import { PropsWithChildren, createContext, useContext, useRef } from "react";
import { loadOmeTiff } from "@hms-dbmi/viv";
import { useChart } from "../../context";
import { createStore } from "zustand";
import { useStoreWithEqualityFn } from "zustand/traditional";
import { observer } from "mobx-react-lite";
import { VivMDVReact } from "../VivMDVReact";
import type { EqFn, Selector, ZustandStore } from "./zustandTypes";

// what about loadOmeZarr, loadBioformatsZarr...
// ... not to mention HTJ2K?
export type OME_TIFF = Awaited<ReturnType<typeof loadOmeTiff>>;

// --- copied straight from Avivator's code::: with notes / changes for MDV ---
import { RENDERING_MODES } from '@hms-dbmi/viv';

const capitalize = string => string.charAt(0).toUpperCase() + string.slice(1);

// typing for generateToggles... not the most useful ones to have, but it's a start.
type BooleanKeys<T> = {
  [K in keyof T]: T[K] extends boolean ? K : never;
}[keyof T];
type TogglesReturnType<T> = {
  [K in BooleanKeys<T> as `toggle${Capitalize<string & K>}`]: () => void;
}
type WithToggles<T> = T & TogglesReturnType<T>;
type SetFunctionType<T> = (fn: (state: T) => T) => void;

function generateToggles<T>(defaults: T, set: SetFunctionType<T>): TogglesReturnType<T> {
  const toggles: any = {};
  Object.entries(defaults).forEach(([k, v]) => {
    if (typeof v === 'boolean') {
      toggles[`toggle${capitalize(k)}`] = () =>
        set(state => ({
          ...state,
          [k]: !state[k]
        }));
    }
  });
  return toggles;
};
export type ChannelsState = {
  channelsVisible: boolean[],
  contrastLimits: [number, number][],
  colors: [number, number, number][],
  domains: [number, number][],
  selections: { z: number, c: number, t: number }[],
  loader: any, //TBD
  image: number,
  ids: string[]
}

export const DEFAUlT_CHANNEL_STATE: ChannelsState = {
  channelsVisible: [] as boolean[],
  contrastLimits: [] as [number, number][],
  colors: [],
  domains: [] as [number, number][],
  selections: [] as { z: number, c: number, t: number }[],
  ids: [],
  // not for serialization... think about this.
  loader: [{ labels: [], shape: [] }],
  image: 0
};
const DEFAUlT_CHANNEL_VALUES = {
  channelsVisible: true,
  contrastLimits: [0, 65535],
  colors: [255, 255, 255],
  domains: [0, 65535],
  selections: { z: 0, c: 0, t: 0 },
  ids: ''
};
const DEFAULT_IMAGE_STATE = {
  lensSelection: 0,
  colormap: '',
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
  onViewportLoad: () => { }
};
type ImageState = typeof DEFAULT_IMAGE_STATE;
const DEFAULT_VIEWER_STATE = {
  isChannelLoading: [] as boolean[],
  isViewerLoading: true,
  pixelValues: [],
  isOffsetsSnackbarOn: false,
  loaderErrorSnackbar: {
    on: false,
    message: null
  },
  isNoImageUrlSnackbarOn: false,
  isVolumeRenderingWarningOn: false,
  useLinkedView: false,
  isControllerOn: true,
  use3d: false,
  useLens: false,
  useColormap: false,
  globalSelection: { z: 0, t: 0 },
  channelOptions: [],
  metadata: null,
  viewState: null,
  source: null as { urlOrFile: string, description: string } | null,
  pyramidResolution: 0
};

// --- following how VivViewerMDV _parseChannels() works ---
export type MdvVivChannelConfig = {
  name: string,
  selections?: { c: number, z: number, t: number }[],
  color?: `#${string}` | [r: number, g: number, b: number],
  visible?: boolean,
  contrastLimits?: [min: number, max: number],
  domains?: [min: number, max: number],
}

export type VivConfig = {
  // todo establish how we actually want this to be...
  channels: MdvVivChannelConfig[],
  // vs...
  image_properties: ChannelsState,
  // vs...
  viewerStore?: Partial<ViewerStore>,
  channelsStore?: Partial<ChannelsStore>,
  imageSettingsStore?: Partial<ImageSettingsStore>,
  //>>>?
  url: string
}

export type ROI = {
  min_x: number,
  min_y: number,
  max_x: number,
  max_y: number,
}

/**
 * In Avivator, there is a zustand store for the channels, image settings, and viewers.
 * While there can be plural 'viewers', it is always of a single image, with shared channels and image settings.
 * 
 * So in order to have appropriate state for multiple MDV charts, we have a context that contains the stores for each.
 * Within this context, we should be able to use much the same API as Avivator.
 */
export type VivContextType = {
  channelsStore: ZustandStore<WithToggles<ChannelsState> & {
    toggleIsOn: (index: number) => void,
    setPropertiesForChannel: (channel: number, newProperties: Partial<typeof DEFAUlT_CHANNEL_VALUES>) => void,
    removeChannel: (channel: number) => void,
    addChannel: (newProperties: Partial<typeof DEFAUlT_CHANNEL_VALUES>) => void
  }>,
  imageSettingsStore: ZustandStore<WithToggles<ImageState>>,
  viewerStore: ZustandStore<WithToggles<typeof DEFAULT_VIEWER_STATE> & {
    setIsChannelLoading: (index: number, val: boolean) => void,
    addIsChannelLoading: (val: boolean) => void,
    removeIsChannelLoading: (index: number) => void
  }> //also has some extra methods... setIsChannelLoading, addIsChannelLoading, removeIsChannelLoading
}
/**
 * This is somewhat MDV-centric now. As of writing, a chart has a property `vivStores`
 * which means different <VivProvider>s (currently rendering in different react roots, but referring to the same chart
 * ie the chart component itself, and a color-change GUI dialog for it)
 * can share reference to the same stores / act as equivalent contexts
 */
export function createVivStores(chart: VivMDVReact) {
  if (chart.vivStores) {
    console.warn('vivStores already exists for this chart');
    return;
  }
  // get any existing values out of the chart's config...
  // currently inactive / experimental state - need to review how props are set,
  // in what way to change these hooks - how much to deviate from Avivator's API
  // (useImage() needs to be changed along with this)
  const channelsStore = createStore(set => ({
    ...DEFAUlT_CHANNEL_STATE,
    // ...chart.config.viv.channelsStore,
    ...generateToggles(DEFAUlT_CHANNEL_VALUES, set),
    toggleIsOn: index =>
      set(state => {
        const channelsVisible = [...state.channelsVisible];
        channelsVisible[index] = !channelsVisible[index];
        return { ...state, channelsVisible };
      }),
    setPropertiesForChannel: (channel, newProperties) =>
      set(state => {
        const entries = Object.entries(newProperties);
        const newState = {};
        entries.forEach(([property, value]) => {
          newState[property] = [...state[property]];
          newState[property][channel] = value;
        });
        return { ...state, ...newState };
      }),
    removeChannel: channel =>
      set(state => {
        const newState = {};
        const channelKeys = Object.keys(DEFAUlT_CHANNEL_VALUES);
        Object.keys(state).forEach(key => {
          if (channelKeys.includes(key)) {
            newState[key] = state[key].filter((_, j) => j !== channel);
          }
        });
        return { ...state, ...newState };
      }),
    addChannel: newProperties =>
      set(state => {
        const entries = Object.entries(newProperties);
        const newState = { ...state };
        entries.forEach(([property, value]) => {
          newState[property] = [...state[property], value];
        });
        Object.entries(DEFAUlT_CHANNEL_VALUES).forEach(([k, v]) => {
          if (newState[k].length < newState[entries[0][0]].length) {
            newState[k] = [...state[k], v];
          }
        });
        return newState;
      })
  }))// satisfies VivContextType['channelsStore'];
  const imageSettingsStore = createStore(set => ({
    ...DEFAULT_IMAGE_STATE,
    // ...chart.config.viv.imageSettingsStore,
    ...generateToggles(DEFAULT_IMAGE_STATE, set)
  })) satisfies VivContextType['imageSettingsStore'];
  const viewerStore = createStore(set => ({
    ...DEFAULT_VIEWER_STATE,
    // ...chart.config.viv.viewerStore,
    ...generateToggles(DEFAULT_VIEWER_STATE, set),
    setIsChannelLoading: (index, val) =>
      set(state => {
        const newIsChannelLoading = [...state.isChannelLoading];
        newIsChannelLoading[index] = val;
        return { ...state, isChannelLoading: newIsChannelLoading };
      }),
    addIsChannelLoading: val =>
      set(state => {
        const newIsChannelLoading = [...state.isChannelLoading, val];
        return { ...state, isChannelLoading: newIsChannelLoading };
      }),
    removeIsChannelLoading: index =>
      set(state => {
        const newIsChannelLoading = [...state.isChannelLoading];
        newIsChannelLoading.splice(index, 1);
        return { ...state, isChannelLoading: newIsChannelLoading };
      })
  })) satisfies VivContextType['viewerStore'];
  // using `as` instead of `satisfies` because of a bug in the type-checking...
  // not sure why 'channelsStore' is misbehaving when other types are ok.
  return { viewerStore, channelsStore, imageSettingsStore } as VivContextType;
}

const VivContext = createContext<VivContextType>(null);
/**
 * This implictly assumes that we are in a context where there is a chart, and that we
 * can access it with `useChart()`. That design decision may be revisited (and it means
 * this detail is deviating somewhat from Avivator).
 * 
 * Separate providers - referring to the same `chart.vivStores` - are used both by the chart
 * itself and by the color-change GUI dialog.
 */
export const VivProvider = observer(({ children }: PropsWithChildren) => {
  const vivChart = useChart() as VivMDVReact; //may want to be less MDV-centric here
  if (!vivChart.vivStores) vivChart.vivStores = createVivStores(vivChart); //<< not very react-y: is this ok?
  return (
    <VivContext.Provider value={vivChart.vivStores}>
      {children}
    </VivContext.Provider>
  )
});

type StoreName = keyof VivContextType;
type ImageSettingsStore = VivContextType['imageSettingsStore'];
type ViewerStore = VivContextType['viewerStore'];
type ChannelsStore = VivContextType['channelsStore'];

// I don't think it's possible to make things like useChannelsStore.setState() work...
// so we have these extra functions...
function useStoreApi<S extends StoreName>(storeName: S): VivContextType[S] {
  const store = useContext(VivContext);
  if (!store) throw 'Viv state hooks must be used within a VivProvider';
  return store[storeName];
}  
export const useChannelsStoreApi = () => useStoreApi('channelsStore');
export const useImageSettingsStoreApi = () => useStoreApi('imageSettingsStore');
export const useViewerStoreApi = () => useStoreApi('viewerStore');

/** should be more-or-less equivalent to equivalent avivator hook - 
 * but there can be multiple viv viewers, so we have context for that.
 */
export function useChannelsStore<U> (selector: Selector<ChannelsStore, U>, equalityFn?: EqFn<U>) {
  const store = useChannelsStoreApi();
  return useStoreWithEqualityFn(store, selector, equalityFn);
}
export function useImageSettingsStore<U>(selector: Selector<ImageSettingsStore, U>, equalityFn?: EqFn<U>) {
  const store = useImageSettingsStoreApi();
  return useStoreWithEqualityFn(store, selector, equalityFn);
}

export function useViewerStore<U>(selector?: Selector<ViewerStore, U>, equalityFn?: EqFn<U>) {
  const store = useViewerStoreApi();
  return useStoreWithEqualityFn(store, selector, equalityFn);
}

export const useLoader = () => {
  const [fullLoader, image] = useChannelsStore(store => [
    store.loader,
    store.image
  ]);
  return Array.isArray(fullLoader[0]) ? fullLoader[image] : fullLoader;
};

export const useMetadata = () => {
  const image = useChannelsStore(store => store.image);
  const metadata = useViewerStore(store => store.metadata);
  return Array.isArray(metadata) ? metadata[image] : metadata;
};

