import * as d3 from "d3";
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
import HistogramWidget, { type HistogramLayer } from "./HistogramWidget";
import { useDebounce } from "use-debounce";

type Range = [number, number];
const HISTOGRAM_BINS = 80;
const HISTOGRAM_WIDTH = 240;
const HISTOGRAM_HEIGHT = 48;

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
                            raster,
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

const clampRange = (range: Range, domain: Range): Range => [
    Math.max(domain[0], Math.min(range[0], domain[1])),
    Math.max(domain[0], Math.min(range[1], domain[1])),
];

const sortRange = ([start, end]: Range): Range =>
    start <= end ? [start, end] : [end, start];

const buildHistogram = (
    values: ArrayLike<number>,
    domain: Range,
    bins: number,
): number[] => {
    const [min, max] = domain;
    if (values.length === 0) {
        return new Array(bins).fill(0);
    }
    const adjustedDomain: Range =
        min === max ? [min, min + 1] : [min, max];
    const histogram = d3
        .bin<number, number>()
        .domain(adjustedDomain)
        .thresholds(bins)(Array.from(values));
    return histogram.map((bin) => bin.length);
};

const ChannelHistogram = ({ index }: { index: number }) => {
    const contrastLimits = useChannelsStore((state) => state.contrastLimits);
    const { colors, domains, raster } = useChannelsStore(({ colors, domains, raster }) => ({
        colors,
        domains,
        raster,
    }));
    const pixelValue = useViewerStore((state) => state.pixelValues[index]);
    const isChannelLoading = useViewerStore((state) => state.isChannelLoading);
    const channelsStore = useChannelsStoreApi();
    const domain = domains[index] ?? ([0, 1] as Range);
    const limits = contrastLimits[index] ?? domain;
    const rasterData = raster[index]?.data;
    const color = colors[index] ?? [37, 99, 235];
    const [liveValue, setLiveValue] = useState<Range | null>(limits);

    useEffect(() => {
        setLiveValue(limits);
    }, [limits]);

    const [debouncedValue] = useDebounce(liveValue, 10, {
        equalityFn: (a, b) => {
            if (!a && !b) return true;
            if (!a || !b) return false;
            return a[0] === b[0] && a[1] === b[1];
        },
    });

    useEffect(() => {
        if (!debouncedValue) return;
        if (debouncedValue.some((value) => Number.isNaN(value))) return;
        if (
            debouncedValue[0] === limits[0] &&
            debouncedValue[1] === limits[1]
        ) {
            return;
        }
        const nextContrastLimits = [...channelsStore.getState().contrastLimits];
        nextContrastLimits[index] = debouncedValue;
        channelsStore.setState({ contrastLimits: nextContrastLimits });
    }, [channelsStore, debouncedValue, index, limits]);

    const layers = useMemo<HistogramLayer[]>(() => {
        const histogramData = buildHistogram(
            rasterData ?? [],
            domain,
            HISTOGRAM_BINS,
        );
        const maxCount = Math.max(1, ...histogramData);
        const highlightData = new Array(HISTOGRAM_BINS).fill(0);
        if (
            Number.isFinite(pixelValue) &&
            domain[0] !== domain[1] &&
            pixelValue >= domain[0] &&
            pixelValue <= domain[1]
        ) {
            const fraction = (pixelValue - domain[0]) / (domain[1] - domain[0]);
            const highlightIndex = Math.min(
                HISTOGRAM_BINS - 1,
                Math.max(0, Math.floor(fraction * HISTOGRAM_BINS)),
            );
            highlightData[highlightIndex] = maxCount;
        }
        return [
            {
                id: `channel-${index}`,
                data: histogramData,
                color: `rgba(${color[0]}, ${color[1]}, ${color[2]}, 0.55)`,
                variant: "bars",
                widthFactor: 0.9,
                inset: 0.05,
                radius: 0.2,
            },
            {
                id: `channel-${index}-pixel`,
                data: highlightData,
                color: "rgba(217, 119, 6, 0.95)",
                variant: "markers",
                widthFactor: 0.35,
                hidden: highlightData.every((value) => value === 0),
            },
        ];
    }, [color, domain, index, pixelValue, rasterData]);

    if (!domains[index] || !contrastLimits[index] || isChannelLoading[index]) {
        return <div className="h-12 px-2 text-sm">Loading...</div>;
    }

    const handleBrushValue = useCallback(
        (value: Range | null) => {
            if (!value) return;
            setLiveValue(sortRange(clampRange(value, domain)));
        },
        [domain],
    );

    const brush = useMemo(
        () => ({
            value: liveValue,
            setValue: handleBrushValue,
            minMax: domain,
        }),
        [domain, handleBrushValue, liveValue],
    );

    return (
        <div className="px-2">
            <HistogramWidget
                layers={layers}
                width={HISTOGRAM_WIDTH}
                height={HISTOGRAM_HEIGHT}
                bins={HISTOGRAM_BINS}
                brush={brush}
            />
        </div>
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
    // not sure I want to be using material.ui... consider adding a widget abstration layer.
    // not jumping right in with using it for all layout etc because I don't want to be tied to it.
    // Starting to use tailwind here.
    // memoizing styles to avoid re-rendering - not sure how to translate to tailwind, not thought about it much yet.
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
                const { domain: domains, contrastLimits, raster } =
                    await getSingleSelectionStats({ loader, selection, use3d });
                const newProps = {
                    domains,
                    contrastLimits,
                    raster,
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
