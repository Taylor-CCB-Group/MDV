import { useCallback, useEffect, useId, useMemo, useState } from "react";
import { shallow } from "zustand/shallow";
import {
    type VivContextType,
    VivProvider,
    useChannelsStore,
    useChannelsStoreApi,
    useLoader,
    useMetadata,
    useViewerStore,
} from "./avivatorish/state";
// some of this is quite bloated - could use dynamic imports for some of the more complex components
import {
    Checkbox,
    FormControl,
    InputLabel,
    MenuItem,
    Select,
    Slider,
} from "@mui/material";
import { PopoverPicker } from "./ColorPicker";
import { getSingleSelectionStats } from "./avivatorish/utils";
import { X } from "lucide-react";
import { useDebounce } from "use-debounce";
import { Histogram, type Range } from "./HistogramComponent";
import { ErrorBoundary } from "react-error-boundary";
import { toJS } from "mobx";
import JsonView from "react18-json-view";
import HistogramWorker from '../../datastore/rawHistogramWorker.ts?worker';
import type { HistogramMessage } from "@/datastore/rawHistogramWorker";

export default function MainVivColorDialog({
    vivStores,
}: { vivStores: VivContextType }) {
    return (
        <VivProvider vivStores={vivStores}>
            <Test />
        </VivProvider>
    );
}

const ChannelChooserMUI = ({ index }: { index: number }) => {
    const channels = useMetadata()?.Pixels.Channels.map(
        (c) => c.Name,
    ) as string[];
    const selections = useChannelsStore(({ selections }) => selections);
    const channelsStore = useChannelsStoreApi();
    const id = useId();
    const name = channels[selections[index].c];

    return (
        <>
            <FormControl fullWidth variant="standard">
                <InputLabel id={id}>{name}</InputLabel>
                <Select
                    labelId={id}
                    label="Channel"
                    onChange={(e) => {
                        const newSelections = [...selections];
                        if (typeof e.target.value !== "number") {
                            return;
                        }
                        newSelections[index].c = e.target.value;
                        channelsStore.setState({ selections: newSelections });
                    }}
                >
                    {channels.map((c, i) => (
                        <MenuItem key={`${i}_${c}`} value={i}>
                            {c}
                        </MenuItem>
                    ))}
                </Select>
            </FormControl>
        </>
    );
};

const ChannelChooser = ({ index }: { index: number }) => {
    const channels = useMetadata()?.Pixels.Channels.map(
        (c) => c.Name,
    ) as string[];
    const { selections, setPropertiesForChannel } = useChannelsStore(
        ({ selections, setPropertiesForChannel }) => ({
            selections,
            setPropertiesForChannel,
        }),
        shallow,
    );
    const loader = useLoader();
    const {
        setIsChannelLoading,
        isChannelLoading,
        removeIsChannelLoading,
        use3d,
    } = useViewerStore(
        ({
            setIsChannelLoading,
            isChannelLoading,
            removeIsChannelLoading,
            use3d,
        }) => ({
            setIsChannelLoading,
            isChannelLoading,
            removeIsChannelLoading,
            use3d,
        }),
        shallow,
    );

    return (
        <>
            <select
                disabled={isChannelLoading[index]}
                value={selections[index].c}
                className="w-full p-1"
                onChange={async (e) => {
                    // see Avivator Controller.jsx onSelectionChange
                    try {
                        const selection = {
                            ...selections[index],
                            c: Number.parseInt(e.target.value),
                        };
                        setIsChannelLoading(index, true);
                        const { domain: domains, contrastLimits, raster } =
                            await getSingleSelectionStats({
                                loader,
                                selection,
                                use3d,
                            });
                        const newProps = {
                            domains,
                            contrastLimits, //, leaving out colors for now - keep existing color
                            raster
                        };
                        setPropertiesForChannel(index, newProps);
                        setIsChannelLoading(index, false);
                        setPropertiesForChannel(index, {
                            selections: selection,
                        });
                    } catch (e) {
                        console.error("failed to load channel");
                        console.error(e);
                        removeIsChannelLoading(index);
                    }
                }}
            >
                {channels.map((c, i) => (
                    <option key={`${i}_${c}`} value={i}>
                        {c}
                    </option>
                ))}
            </select>
        </>
    );
};

const BrightnessContrast = ({ index }: { index: number }) => {
    const { contrast, brightness } = useChannelsStore(
        ({ contrast, brightness }) => ({ contrast, brightness }),
    );
    const channelsStore = useChannelsStoreApi();
    const isChannelLoading = useViewerStore((state) => state.isChannelLoading);
    return (
        <div className="col-span-5 pl-2 pr-2 flex items-center gap-4">
            contrast:
            <Slider
                size="small"
                disabled={isChannelLoading[index]}
                value={contrast[index]}
                min={0.1}
                max={0.9}
                step={0.01}
                onChange={(_, v) => {
                    contrast[index] = v as number;
                    const newContrast = [...contrast];
                    channelsStore.setState({ contrast: newContrast });
                }}
                onClick={(e) => {
                    // reset value on double-click
                    if (e.detail === 2) {
                        contrast[index] = 0.5;
                        const newContrast = [...contrast];
                        channelsStore.setState({ contrast: newContrast });
                    }
                }}
            />
            brightness:
            <Slider
                size="small"
                disabled={isChannelLoading[index]}
                value={brightness[index]}
                min={0.01}
                max={0.99}
                step={0.01}
                onChange={(_, v) => {
                    brightness[index] = v as number;
                    const newbrightness = [...brightness];
                    channelsStore.setState({ brightness: newbrightness });
                }}
                onClick={(e) => {
                    // reset value on double-click
                    if (e.detail === 2) {
                        brightness[index] = 0.5;
                        const newbrightness = [...brightness];
                        channelsStore.setState({ brightness: newbrightness });
                    }
                }}
            />
        </div>
    );
};

/** 
 * @param index - index of the channel in the channels store - NOT the index of the channel in the metadata
 * @returns the name of the channel
 */
const useChannelName = (index: number) => {
    const selections = useChannelsStore(({ selections }) => selections, shallow);
    const { c } = selections[index];
    return useMetadata()?.Pixels.Channels[c].Name;
}

const ChannelHistogram = ({ index }: { index: number }) => {
    // return <div />
    //using shallow as per Avivator *prevents* re-rendering which should be happening << review
    const limits = useChannelsStore(({ contrastLimits }) => contrastLimits);
    // should be 'rasters' really
    const { domains, raster } = useChannelsStore(({ domains, raster }) => ({ domains, raster }));
    const { pixelValues } = useViewerStore(({ pixelValues }) => ({ pixelValues }));
    const pixelValue = pixelValues[index];
    const domain = domains[index];
    const rasterData = raster[index]?.data || [0]; //! revisit this sometime

    const [histogramData, setHistogramData] = useState([] as number[]);

    // review this / use d3 scales, see if we can simplify state
    const [min, max] = domain;
    const scaleValue = useCallback((v: number) => (v - min) / (max - min), [min, max]);
    const { limit, normalisedLow, normalisedHigh } = useMemo(() => {
        const limit = limits[index];
        // this can probably be integrated into the histogram component
        const normalisedLow = scaleValue(limit[0]);
        const normalisedHigh = scaleValue(limit[1]);
        // const normalisedPixelValue = scaleValue(pixelValue);
        return { limit, normalisedLow, normalisedHigh };
    }, [limits, index, scaleValue]);
    const channelsStore = useChannelsStoreApi();


    // todo review whether this state can be cleaned up / simplified
    const [liveValue, setLiveValue] = useState<Range>([0, 0]);
    const debounceTime = 10; //how low can we go?
    const [debouncedValue] = useDebounce(liveValue, debounceTime);
    useEffect(() => {
        // this feels glitchy, need to iron out some issues
        // maybe use useRef to store the value?
        if (!debouncedValue) return;
        // wtf, why is this needed?
        // maybe related to setValue(null), which we do sometimes when clearing the brush
        // but it seems to be happening during other edits too?
        if (debouncedValue.some(Number.isNaN)) return;
        if (debouncedValue[0] === 0 && debouncedValue[1] === 0) return;
        // debouncing won't help if we still have a dependency on limits
        const limits = channelsStore.getState().contrastLimits;
        if (debouncedValue[0] === limits[index][0] && debouncedValue[1] === limits[index][1]) return;
        limits[index] = debouncedValue;
        const contrastLimits = [...limits];
        // we have a weird mix of mobx and zustand here - don't think that's a particular source of bug,
        // but it smells quite a bit
        channelsStore.setState({ contrastLimits });
    }, [debouncedValue, index, channelsStore]);
    const channelName = useChannelName(index);
    const isChannelLoading = useViewerStore((state) => state.isChannelLoading)[index];
    if (isChannelLoading) return <div>Loading...</div>;
    return (
        <div className="p-4">
            {/* {isChannelLoading ?  */}
            <Histogram
                value={limit}
                step={0.001} //!todo make this dynamic
                data={rasterData}
                // todo add value label
                lowFraction={normalisedLow} // component should calculate these
                highFraction={normalisedHigh}
                setValue={v => {
                    // we want a different behaviour for clearing the brush here vs selection dialog...
                    if (v) setLiveValue(v);
                }}
                domain={domain}
                bins={150}
                histoHeight={50}
                highlightValue={pixelValue}
                //todo make scales adapt or be configurable
                //log scale breaks with 0 values
                xScaleType="linear"
                yScaleType="symlog"
                name={channelName}
            />
            {/* : "Loading..."} */}
        </div>
    )
}

const ChannelController = ({ index }: { index: number }) => {
    const { colors, channelsVisible, removeChannel } =
        useChannelsStore(
            ({ colors, channelsVisible, removeChannel }) => ({
                colors,
                channelsVisible,
                removeChannel,
            }),
        );
    const isChannelLoading = useViewerStore((state) => state.isChannelLoading);
    const metadata = useMetadata();
    const channelsStore = useChannelsStoreApi();
    if (!metadata) throw "no metadata"; //TODO type metadata
    const channelVisible = channelsVisible[index];
    const color = colors[index];
    const gridStyle = useMemo(
        () => ({
            gridTemplateColumns: "0.4fr 0.1fr 0.1fr 1fr 0.1fr",
        }),
        [],
    );
    return (
        <>
            <div
                className="grid justify-start items-center p-1"
                style={gridStyle}
            >
                <ChannelChooser index={index} />
                <Checkbox
                    checked={channelVisible}
                    disabled={isChannelLoading[index]}
                    onChange={() => {
                        channelsVisible[index] = !channelVisible;
                        const visible = [...channelsVisible];
                        channelsStore.setState({ channelsVisible: visible });
                    }}
                />
                <PopoverPicker
                    color={color}
                    onChange={(c) => {
                        colors[index] = c;
                        const newColors = [...colors];
                        channelsStore.setState({ colors: newColors });
                    }}
                />
                <ErrorBoundary fallback={<div>Histogram error</div>}>
                    <ChannelHistogram index={index} />
                </ErrorBoundary>
                <button
                    type="button"
                    className="pl-4"
                    onClick={() => {
                        removeChannel(index);
                    }}
                >
                    <X />
                </button>
                <BrightnessContrast index={index} />
            </div>
        </>
    );
};

const AddChannel = () => {
    const loader = useLoader();
    // const { labels } = loader[0];
    // const channelsStore = useChannelsStoreApi();
    const { selections, setPropertiesForChannel } = useChannelsStore(
        ({ selections, setPropertiesForChannel }) => ({
            selections,
            setPropertiesForChannel,
        }),
    );
    const canAddChannel = selections.length < 6;
    const addChannel = useChannelsStore((state) => state.addChannel);
    const { use3d, setIsChannelLoading } = useViewerStore(
        ({ use3d, setIsChannelLoading }) => ({ use3d, setIsChannelLoading }),
        shallow,
    );
    return (
        <button
            type="button"
            className="p-2 rounded-lg bg-slate-400"
            disabled={!canAddChannel}
            onClick={async () => {
                // would be nice to have less repition of this code here and in ChannelController
                const index = selections.length;
                setIsChannelLoading(index, true);
                const selection = { c: 0, z: 0, t: 0 };
                addChannel({
                    selections: selection,
                    ids: String(Math.random()),
                    channelsVisible: true,
                });
                const { domain: domains, contrastLimits } =
                    await getSingleSelectionStats({ loader, selection, use3d });
                const newProps = {
                    domains,
                    contrastLimits,
                    //raster
                };
                setPropertiesForChannel(index, newProps);
                setIsChannelLoading(index, false);
            }}
        >
            Add channel
        </button>
    );
};
function ChannelError() {
    //this can be simpler with non-Api version
    const store = useChannelsStoreApi();
    const state = useMemo(() => {
        return toJS(store.getState());
    }, [store]);
    return (
        <JsonView src={state} />
    );
}
export const Test = () => {
    // ids are strings like "0.1234", and somehow there are too many of them
    // when image was saved with removed channels - e.g. default 4 ids appear when
    // image is saved with 2 channels
    // seems like an inconsistency in the store
    const ids = useChannelsStore(({ ids }) => ids);
    const selections = useChannelsStore(({ selections }) => selections);
    if (selections.length !== ids.length) {
        console.error(
            "selections and ids do not match",
            selections,
            ids,
        );
    }
    return (
        <div className="bg-[hsl(var(--input))] w-full">
            {ids.map((id, i) => (
                <ErrorBoundary fallback={<ChannelError />} key={id}>
                    <ChannelController key={id} index={i} />
                </ErrorBoundary>
            ))}
            <AddChannel />
        </div>
    );
};
