import { useCallback, useEffect, useId, useMemo, useState } from "react";
import { shallow } from "zustand/shallow";
import {
    type VivContextType,
    VivProvider,
    useChannelsStore,
    useChannelsStoreApi,
    useImageSettingsStoreApi,
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

const ChannelHistogram = ({ index }: { index: number }) => {
    // return <div />
    const limits = useChannelsStore(({ contrastLimits }) => contrastLimits); //using shallow as per Avivator *prevents* re-rendering which should be happening
    // should be 'rasters' really
    const { domains, raster } = useChannelsStore(({ domains, raster }) => ({ domains, raster }));
    const { pixelValues } = useViewerStore(({ pixelValues }) => ({ pixelValues }));
    const pixelValue = pixelValues[index];
    const domain = domains[index];
    const rasterData = raster[index]?.data || [0]; //! revisit this sometime
    
    const [histogramData, setHistogramData] = useState([] as number[]);

    const [min, max] = domain;
    const scaleValue = useCallback((v: number) => (v - min) / (max - min), [min, max]);
    const normalisedPixelValue = scaleValue(pixelValue);
    const limit = useMemo(() => limits[index], [limits, index]);
    const normalisedLow = scaleValue(limit[0]);
    const normalisedHigh = scaleValue(limit[1]);
    const channelsStore = useChannelsStoreApi();


    // this is never being called?
    const queryHistogram = useCallback(async () => {
        // todo nicer worker syntax?
        const worker = new Worker(
            new URL("../../datastore/rawHistogramWorker.ts", import.meta.url),
        );
        worker.onmessage = (event) => {
            setHistogramData(event.data);
            worker.terminate();
        };
        //! this is a lie - need to actually check!
        const isInt32 = false;
        const originalData = rasterData;
        const data = new SharedArrayBuffer(originalData.length * 4);
        new Float32Array(data).set(originalData);
        worker.postMessage({
            data,
            min,
            max,
            bins: 100,
            isFloat: isInt32,
        });
    }, [min, max, rasterData]);
    const [liveValue, setLiveValue] = useState<Range>([0, 0]);
    const [debouncedValue] = useDebounce(liveValue, 10);
    useEffect(() => {
        // this feels glitchy, need to iron out some issues
        if (!debouncedValue) return;
        if (debouncedValue[0] === 0 && debouncedValue[1] === 0) return;
        // debouncing won't help if we still have a dependency on limits
        const limits = channelsStore.getState().contrastLimits;
        if (debouncedValue[0] === limits[index][0] && debouncedValue[1] === limits[index][1]) return;
        limits[index] = debouncedValue;
        const contrastLimits = [...limits];
        channelsStore.setState({ contrastLimits });
    }, [debouncedValue, index, channelsStore]);

    return (
        <Histogram 
            value={limit}
            step={0.001} //todo make this dynamic
            histogram={histogramData}
            // todo add indicator for highightValue={pixelValue}
            lowFraction={normalisedLow} // component should calculate these
            highFraction={normalisedHigh}
            queryHistogram={queryHistogram}
            setValue={setLiveValue}
            minMax={domain}
            histoWidth={100}
            histoHeight={50}
        />
    )
}

const ChannelController = ({ index }: { index: number }) => {
    const limits = useChannelsStore(({ contrastLimits }) => contrastLimits); //using shallow as per Avivator *prevents* re-rendering which should be happening
    const { colors, domains, channelsVisible, removeChannel } =
        useChannelsStore(
            ({ colors, domains, channelsVisible, removeChannel }) => ({
                colors,
                domains,
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
    const colorString = useMemo(() => `rgb(${color[0]}, ${color[1]}, ${color[2]})`, [color[0], color[1], color[2]]);
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
                <ChannelHistogram index={index} />
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

export const Test = () => {
    const ids = useChannelsStore(({ ids }) => ids);
    return (
        <div className="bg-[hsl(var(--input))] w-full">
            {ids.map((id, i) => (
                <ChannelController key={id} index={i} />
            ))}
            <AddChannel />
        </div>
    );
};
