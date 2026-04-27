import { resolveAutoHistogramXScaleFromValues, resolveAutoHistogramYScale } from "@/lib/utils";
import { isArray } from "@/lib/utils";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import HighlightOffIcon from "@mui/icons-material/HighlightOff";
// some of this is quite bloated - could use dynamic imports for some of the more complex components
import {
    Accordion,
    AccordionDetails,
    AccordionSummary,
    Checkbox,
    FormControl,
    IconButton,
    InputLabel,
    MenuItem,
    Select,
    Slider,
} from "@mui/material";
import * as d3 from "d3";
import { useCallback, useEffect, useId, useMemo, useRef, useState } from "react";
import { useDebounce } from "use-debounce";
import { shallow } from "zustand/shallow";
import { useTheme } from "../hooks";
import { PopoverPicker } from "./ColorPicker";
import HistogramWidget, { type HistogramLayer, type HistogramScaleType } from "./HistogramWidget";
import {
    type VivContextType,
    VivProvider,
    useChannelsStore,
    useChannelsStoreApi,
    useLoader,
    useMetadata,
    useViewerStore,
} from "./avivatorish/state";
import { getSingleSelectionStats } from "./avivatorish/utils";

type Range = [number, number];
type ScaleMode = "auto" | HistogramScaleType;
const HISTOGRAM_BINS = 200;
const HISTOGRAM_WIDTH = 240;
const HISTOGRAM_HEIGHT = 78;
const stopAccordionToggle = (event: React.SyntheticEvent) => {
    event.stopPropagation();
};
const toggleScaleMode = (mode: ScaleMode, resolvedMode: HistogramScaleType): ScaleMode => {
    if (mode === "auto") {
        return resolvedMode === "log" ? "linear" : "log";
    }
    return mode === "log" ? "linear" : "log";
};

export default function MainVivColorDialog({ vivStores }: { vivStores: VivContextType }) {
    return (
        <VivProvider vivStores={vivStores}>
            <Test />
        </VivProvider>
    );
}

const ChannelChooserMUI = ({ index }: { index: number }) => {
    const channels = useMetadata()?.Pixels.Channels.map((c) => c.Name) as string[];
    const selections = useChannelsStore(({ selections }) => selections);
    const channelsStore = useChannelsStoreApi();
    const id = useId();
    const name = channels[selections[index].c ?? 0];

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
    const channels = useMetadata()?.Pixels.Channels.map((c) => c.Name) as string[];
    const { selections, setPropertiesForChannel } = useChannelsStore(
        ({ selections, setPropertiesForChannel }) => ({
            selections,
            setPropertiesForChannel,
        }),
        shallow,
    );
    const loader = useLoader();
    const { setIsChannelLoading, isChannelLoading, removeIsChannelLoading, use3d } = useViewerStore(
        ({ setIsChannelLoading, isChannelLoading, removeIsChannelLoading, use3d }) => ({
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
                value={selections[index].c ?? 0}
                className="w-full rounded-md border border-[hsl(var(--border))] bg-[hsl(var(--background))] px-2 py-1.5 text-sm text-[hsl(var(--foreground))] shadow-sm outline-none transition focus:border-[hsl(var(--ring))] disabled:cursor-wait disabled:opacity-60"
                onChange={async (e) => {
                    // see Avivator Controller.jsx onSelectionChange
                    try {
                        const selection = {
                            ...selections[index],
                            c: Number.parseInt(e.target.value),
                        };
                        setIsChannelLoading(index, true);
                        const {
                            domain: domains,
                            contrastLimits,
                            raster,
                        } = await getSingleSelectionStats({
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

const sortRange = ([start, end]: Range): Range => (start <= end ? [start, end] : [end, start]);

const rangesEqual = (a: Range, b: Range) => a[0] === b[0] && a[1] === b[1];

const buildHistogram = (
    values: ArrayLike<number>,
    domain: Range,
    bins: number,
    xScaleType: HistogramScaleType,
): { counts: number[]; edges: number[] } => {
    const [min, max] = domain;
    if (values.length === 0) {
        return {
            counts: new Array(bins).fill(0),
            edges: Array.from({ length: bins + 1 }, (_, index) => min + ((max - min) * index) / bins),
        };
    }
    const adjustedDomain: Range = min === max ? [min, min + 1] : [min, max];
    const baseHistogram = d3.bin<number, number>().domain(adjustedDomain);
    const histogram =
        xScaleType === "log"
            ? baseHistogram.thresholds(
                  Array.from({ length: bins - 1 }, (_, index) => {
                      const t = (index + 1) / bins;
                      return d3.scaleSymlog().domain(adjustedDomain).range([0, 1]).invert(t);
                  }),
              )(Array.from(values))
            : baseHistogram.thresholds(bins)(Array.from(values));
    const counts = histogram.map((bin) => bin.length);
    const edges = [histogram[0]?.x0 ?? adjustedDomain[0], ...histogram.map((bin) => bin.x1 ?? adjustedDomain[1])];
    return { counts, edges };
};

const ChannelHistogram = ({ index }: { index: number }) => {
    const contrastLimits = useChannelsStore((state) => state.contrastLimits);
    const color = useChannelsStore((state) => state.colors[index] ?? [37, 99, 235]);
    const currentDomain = useChannelsStore((state) => state.domains[index]);
    const domain = currentDomain ?? ([0, 1] as Range);
    const rasterData = useChannelsStore((state) => state.raster[index]?.data);
    const pixelValue = useViewerStore((state) => state.pixelValues[index]);
    const isChannelLoading = useViewerStore((state) => state.isChannelLoading);
    const channelsStore = useChannelsStoreApi();
    const dark = useTheme() === "dark";
    const limits = contrastLimits[index] ?? domain;
    const [liveValue, setLiveValue] = useState<Range | null>(limits);
    const [xScaleMode, setXScaleMode] = useState<ScaleMode>("auto");
    const [yScaleMode, setYScaleMode] = useState<ScaleMode>("auto");

    useEffect(() => {
        setLiveValue(rangesEqual(limits, domain) ? null : limits);
    }, [domain, limits]);

    const [debouncedValue] = useDebounce(liveValue, 10, {
        equalityFn: (a, b) => {
            if (!a && !b) return true;
            if (!a || !b) return false;
            return a[0] === b[0] && a[1] === b[1];
        },
    });

    useEffect(() => {
        const nextValue = debouncedValue ?? domain;
        if (nextValue.some((value) => Number.isNaN(value))) return;
        if (rangesEqual(nextValue, limits)) {
            return;
        }
        const nextContrastLimits = [...channelsStore.getState().contrastLimits];
        nextContrastLimits[index] = nextValue;
        channelsStore.setState({ contrastLimits: nextContrastLimits });
    }, [channelsStore, debouncedValue, domain, index, limits]);

    const resolvedXScale = useMemo(
        () => (xScaleMode === "auto" ? resolveAutoHistogramXScaleFromValues(domain, rasterData) : xScaleMode),
        [domain, rasterData, xScaleMode],
    );
    const histogram = useMemo(
        () => buildHistogram(rasterData ?? [], domain, HISTOGRAM_BINS, resolvedXScale),
        [domain, rasterData, resolvedXScale],
    );
    const resolvedYScale = useMemo(
        () => (yScaleMode === "auto" ? resolveAutoHistogramYScale(histogram.counts) : yScaleMode),
        [histogram.counts, yScaleMode],
    );

    const layers = useMemo<HistogramLayer[]>(() => {
        return [
            {
                id: `channel-${index}`,
                data: histogram.counts,
                color: `rgba(${color[0]}, ${color[1]}, ${color[2]}, 0.55)`,
                variant: "line",
            },
        ];
    }, [color, histogram.counts, index]);

    const isHistogramLoading = !currentDomain || !contrastLimits[index] || isChannelLoading[index];

    const handleBrushValue = useCallback(
        (value: Range | null) => {
            if (!value) {
                setLiveValue(null);
                return;
            }
            const nextValue = sortRange(clampRange(value, domain));
            setLiveValue(rangesEqual(nextValue, domain) ? null : nextValue);
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
        <div className="flex h-full min-h-[94px] rounded-md border border-[hsl(var(--border))] bg-[hsl(var(--background))] p-1.5">
            {isHistogramLoading ? (
                <div className="flex h-[68px] w-full items-center rounded-md border border-dashed border-[hsl(var(--border))] bg-[hsl(var(--muted))] px-2 text-[10px] font-medium uppercase tracking-[0.14em] text-[hsl(var(--muted-foreground))]">
                    Loading histogram
                </div>
            ) : (
                <HistogramWidget
                    layers={layers}
                    width={HISTOGRAM_WIDTH}
                    height={HISTOGRAM_HEIGHT}
                    bins={HISTOGRAM_BINS}
                    binEdges={histogram.edges}
                    xScaleType={resolvedXScale}
                    yScaleType={resolvedYScale}
                    brush={brush}
                    scaleControls={{
                        xLabel: xScaleMode === "auto" ? resolvedXScale : xScaleMode,
                        yLabel: yScaleMode === "auto" ? resolvedYScale : yScaleMode,
                        onToggleX: () => setXScaleMode((mode) => toggleScaleMode(mode, resolvedXScale)),
                        onToggleY: () => setYScaleMode((mode) => toggleScaleMode(mode, resolvedYScale)),
                    }}
                    markers={[
                        {
                            id: `pixel-value-${index}`,
                            value: pixelValue,
                            color: dark ? "rgba(255, 255, 255, 0.8)" : "rgba(15, 23, 42, 0.8)",
                            hidden: !Number.isFinite(pixelValue) || pixelValue < domain[0] || pixelValue > domain[1],
                        },
                    ]}
                />
            )}
        </div>
    );
};

const BrightnessContrast = ({ index }: { index: number }) => {
    const { contrast, brightness } = useChannelsStore(({ contrast, brightness }) => ({ contrast, brightness }));
    const channelsStore = useChannelsStoreApi();
    const isChannelLoading = useViewerStore((state) => state.isChannelLoading);
    return (
        <div className="grid gap-x-4 gap-y-2 sm:grid-cols-2">
            <div className="flex items-center gap-2">
                <span className="w-16 shrink-0 text-[10px] font-medium uppercase tracking-[0.14em] text-[hsl(var(--muted-foreground))]">
                    Contrast
                </span>
                <Slider
                    size="small"
                    disabled={isChannelLoading[index]}
                    value={contrast[index]}
                    min={0.1}
                    max={0.9}
                    step={0.01}
                    onChange={(_, v) => {
                        if (isArray(v)) return;
                        const newContrast = [...contrast];
                        newContrast[index] = v;
                        channelsStore.setState({ contrast: newContrast });
                    }}
                    onClick={(e) => {
                        if (e.detail === 2) {
                            const newContrast = [...contrast];
                            newContrast[index] = 0.5;
                            channelsStore.setState({ contrast: newContrast });
                        }
                    }}
                />
            </div>
            <div className="flex items-center gap-2">
                <span className="w-16 shrink-0 text-[10px] font-medium uppercase tracking-[0.14em] text-[hsl(var(--muted-foreground))]">
                    Brightness
                </span>
                <Slider
                    size="small"
                    disabled={isChannelLoading[index]}
                    value={brightness[index]}
                    min={0.01}
                    max={0.99}
                    step={0.01}
                    onChange={(_, v) => {
                        if (isArray(v)) return;
                        const newBrightness = [...brightness];
                        newBrightness[index] = v;
                        channelsStore.setState({ brightness: newBrightness });
                    }}
                    onClick={(e) => {
                        if (e.detail === 2) {
                            const newBrightness = [...brightness];
                            newBrightness[index] = 0.5;
                            channelsStore.setState({ brightness: newBrightness });
                        }
                    }}
                />
            </div>
        </div>
    );
};

const ChannelController = ({ index }: { index: number }) => {
    const color = useChannelsStore((state) => state.colors[index]);
    const channelVisible = useChannelsStore((state) => state.channelsVisible[index]);
    const removeChannel = useChannelsStore((state) => state.removeChannel);
    const isChannelLoading = useViewerStore((state) => state.isChannelLoading);
    const metadata = useMetadata();
    const channelsStore = useChannelsStoreApi();
    const [isHovered, setIsHovered] = useState(false);

    if (!metadata) throw "no metadata"; //TODO type metadata
    if (!color) return null;

    const hasPendingLoads = isChannelLoading.some(Boolean);

    return (
        <Accordion
            disableGutters
            defaultExpanded={false}
            elevation={0}
            onMouseEnter={() => setIsHovered(true)}
            onMouseLeave={() => setIsHovered(false)}
            sx={{
                position: "relative",
                border: "1px solid hsl(var(--border))",
                borderRadius: "8px",
                backgroundColor: "hsl(var(--background))",
                boxShadow: "0 1px 2px hsl(var(--foreground) / 0.04)",
                "&:before": { display: "none" },
            }}
        >
            <AccordionSummary
                expandIcon={<ExpandMoreIcon fontSize="small" />}
                sx={{
                    px: 1.25,
                    py: 0.75,
                    minHeight: "unset",
                    "& .MuiAccordionSummary-content": {
                        my: 0,
                        minWidth: 0,
                    },
                    "& .MuiAccordionSummary-expandIconWrapper": {
                        color: "hsl(var(--muted-foreground))",
                    },
                }}
            >
                <div
                    className="grid w-full min-w-0 grid-cols-[minmax(11rem,0.95fr)_minmax(0,1.15fr)] gap-2"
                    onClick={stopAccordionToggle}
                    onFocus={stopAccordionToggle}
                >
                    <div className="flex min-h-[86px] min-w-0 flex-col justify-between rounded-md border border-[hsl(var(--border))] bg-[hsl(var(--background))] px-2 py-1.5">
                        <div className="flex items-center justify-between gap-2 text-[10px] font-medium uppercase tracking-[0.14em] text-[hsl(var(--muted-foreground))]">
                            <div className="flex items-center gap-2">
                                <span
                                    className="h-2.5 w-2.5 rounded-full border border-black/10"
                                    style={{
                                        backgroundColor: `rgb(${color[0]}, ${color[1]}, ${color[2]})`,
                                    }}
                                />
                                <span>Channel</span>
                            </div>
                            <div className="flex items-center gap-1.5">
                                <Checkbox
                                    checked={channelVisible}
                                    disabled={isChannelLoading[index]}
                                    size="small"
                                    sx={{ padding: "2px" }}
                                    onClick={stopAccordionToggle}
                                    onFocus={stopAccordionToggle}
                                    onChange={() => {
                                        const channelsVisible = channelsStore.getState().channelsVisible;
                                        const visible = [...channelsVisible];
                                        visible[index] = !channelVisible;
                                        channelsStore.setState({ channelsVisible: visible });
                                    }}
                                />
                                <div
                                    className="h-7 min-w-9 overflow-hidden rounded border border-[hsl(var(--border))] bg-[hsl(var(--background))] shadow-sm"
                                    onClick={stopAccordionToggle}
                                    onFocus={stopAccordionToggle}
                                >
                                    <PopoverPicker
                                        color={color}
                                        onChange={(c) => {
                                            const colors = channelsStore.getState().colors;
                                            const newColors = [...colors];
                                            newColors[index] = c;
                                            channelsStore.setState({ colors: newColors });
                                        }}
                                    />
                                </div>
                            </div>
                        </div>
                        <div onClick={stopAccordionToggle} onFocus={stopAccordionToggle}>
                            <ChannelChooser index={index} />
                        </div>
                    </div>
                    <div className="min-w-0 self-stretch">
                        <ChannelHistogram index={index} />
                    </div>
                </div>
                <IconButton
                    disabled={hasPendingLoads}
                    onClick={(event) => {
                        event.stopPropagation();
                        removeChannel(index);
                    }}
                    size="small"
                    aria-label="remove channel"
                    sx={{
                        position: "absolute",
                        right: "-14px",
                        top: "-14px",
                        color: "hsl(var(--muted-foreground))",
                        opacity: isHovered ? 1 : 0,
                        transition: "opacity 0.3s",
                        "&:focus": {
                            opacity: 1,
                        },
                    }}
                >
                    <HighlightOffIcon fontSize="small" />
                </IconButton>
            </AccordionSummary>
            <AccordionDetails
                sx={{
                    px: 1.25,
                    pb: 1.25,
                    pt: 0,
                }}
            >
                <div className="rounded-md border border-[hsl(var(--border))] bg-[hsl(var(--muted) / 0.3)] px-2 py-2">
                    <div className="mb-2 text-[10px] font-medium uppercase tracking-[0.14em] text-[hsl(var(--muted-foreground))]">
                        Tone controls
                    </div>
                    <BrightnessContrast index={index} />
                </div>
            </AccordionDetails>
        </Accordion>
    );
};

const AddChannel = () => {
    const loader = useLoader();
    const isAddingChannelRef = useRef(false);
    const [isAddingChannel, setIsAddingChannel] = useState(false);
    // const { labels } = loader[0];
    // const channelsStore = useChannelsStoreApi();
    const { selections, setPropertiesForChannel } = useChannelsStore(({ selections, setPropertiesForChannel }) => ({
        selections,
        setPropertiesForChannel,
    }));
    const canAddChannel = selections.length < 6;
    const { addChannel, removeChannel } = useChannelsStore((state) => ({
        addChannel: state.addChannel,
        removeChannel: state.removeChannel,
    }));
    const { use3d, setIsChannelLoading, removeIsChannelLoading } = useViewerStore(
        ({ use3d, setIsChannelLoading, removeIsChannelLoading }) => ({
            use3d,
            setIsChannelLoading,
            removeIsChannelLoading,
        }),
        shallow,
    );
    return (
        <button
            type="button"
            className="rounded-md border border-[hsl(var(--border))] bg-[hsl(var(--background))] px-3 py-1.5 text-sm font-medium text-[hsl(var(--foreground))] shadow-sm transition hover:bg-[hsl(var(--muted)/0.45)] disabled:cursor-not-allowed disabled:opacity-40"
            disabled={!canAddChannel || isAddingChannel}
            onClick={async () => {
                if (isAddingChannelRef.current) {
                    return;
                }
                isAddingChannelRef.current = true;
                setIsAddingChannel(true);
                // would be nice to have less repition of this code here and in ChannelController
                const index = selections.length;
                let didAddChannel = false;
                try {
                    setIsChannelLoading(index, true);
                    const baseSelection = selections[0] ?? {};
                    const selection = {
                        ...baseSelection,
                        c: baseSelection.c ?? 0,
                    };
                    addChannel({
                        selections: selection,
                        ids: String(Math.random()),
                        channelsVisible: true,
                    });
                    didAddChannel = true;
                    const {
                        domain: domains,
                        contrastLimits,
                        raster,
                    } = await getSingleSelectionStats({ loader, selection, use3d });
                    const newProps = {
                        domains,
                        contrastLimits,
                        raster,
                    };
                    setPropertiesForChannel(index, newProps);
                    setIsChannelLoading(index, false);
                } catch (error) {
                    console.error("failed to add channel");
                    console.error(error);
                    if (didAddChannel) {
                        removeChannel(index);
                        removeIsChannelLoading(index);
                        return;
                    }
                    setIsChannelLoading(index, false);
                } finally {
                    isAddingChannelRef.current = false;
                    setIsAddingChannel(false);
                }
            }}
        >
            Add channel
        </button>
    );
};

export const Test = () => {
    const ids = useChannelsStore(({ ids }) => ids);
    return (
        <div className="w-full space-y-2 bg-[hsl(var(--background))] p-2">
            {ids.map((id, i) => (
                <ChannelController key={id} index={i} />
            ))}
            <AddChannel />
        </div>
    );
};
