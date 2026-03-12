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
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import HighlightOffIcon from "@mui/icons-material/HighlightOff";
import { PopoverPicker } from "./ColorPicker";
import { getSingleSelectionStats } from "./avivatorish/utils";
import HistogramWidget, {
    type HistogramLayer,
    type HistogramScaleType,
} from "./HistogramWidget";
import { useDebounce } from "use-debounce";

type Range = [number, number];
type ScaleMode = "auto" | HistogramScaleType;
const HISTOGRAM_BINS = 80;
const HISTOGRAM_WIDTH = 240;
const HISTOGRAM_HEIGHT = 78;
const toggleScaleMode = (
    mode: ScaleMode,
    resolvedMode: HistogramScaleType,
): ScaleMode => {
    if (mode === "auto") {
        return resolvedMode === "log" ? "linear" : "log";
    }
    return mode === "log" ? "linear" : "log";
};

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
                className="w-full rounded-md border border-[hsl(var(--border))] bg-[hsl(var(--background))] px-2 py-1.5 text-sm text-[hsl(var(--foreground))] shadow-sm outline-none transition focus:border-[hsl(var(--ring))] disabled:cursor-wait disabled:opacity-60"
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

const resolveAutoXScale = (domain: Range): HistogramScaleType => {
    const [min, max] = domain;
    if (!Number.isFinite(min) || !Number.isFinite(max) || min === max) {
        return "linear";
    }
    const shiftedMin = Math.min(Math.abs(min), Math.abs(max)) < 1e-9
        ? 1e-9
        : Math.max(1e-9, Math.min(Math.abs(min), Math.abs(max)));
    const shiftedMax = Math.max(Math.abs(min), Math.abs(max), shiftedMin);
    return shiftedMax / shiftedMin > 500 ? "log" : "linear";
};

const resolveAutoYScale = (histogram: number[]): HistogramScaleType => {
    const nonZero = histogram.filter((value) => value > 0);
    if (nonZero.length < 2) return "linear";
    const min = Math.min(...nonZero);
    const max = Math.max(...nonZero);
    const mean = nonZero.reduce((sum, value) => sum + value, 0) / nonZero.length;
    return max / min > 50 || max / Math.max(1, mean) > 10 ? "log" : "linear";
};

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
            edges: Array.from({ length: bins + 1 }, (_, index) =>
                min + ((max - min) * index) / bins,
            ),
        };
    }
    const adjustedDomain: Range =
        min === max ? [min, min + 1] : [min, max];
    const baseHistogram = d3.bin<number, number>().domain(adjustedDomain);
    const histogram =
        xScaleType === "log"
            ? baseHistogram.thresholds(
                  Array.from({ length: bins - 1 }, (_, index) => {
                      const t = (index + 1) / bins;
                      return d3
                          .scaleSymlog()
                          .domain(adjustedDomain)
                          .range([0, 1])
                          .invert(t);
                  }),
              )(Array.from(values))
            : baseHistogram.thresholds(bins)(Array.from(values));
    const counts = histogram.map((bin) => bin.length);
    const edges = [
        histogram[0]?.x0 ?? adjustedDomain[0],
        ...histogram.map((bin) => bin.x1 ?? adjustedDomain[1]),
    ];
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
    const limits = contrastLimits[index] ?? domain;
    const [liveValue, setLiveValue] = useState<Range | null>(limits);
    const [xScaleMode, setXScaleMode] = useState<ScaleMode>("auto");
    const [yScaleMode, setYScaleMode] = useState<ScaleMode>("auto");

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

    const resolvedXScale = useMemo(
        () => (xScaleMode === "auto" ? resolveAutoXScale(domain) : xScaleMode),
        [domain, xScaleMode],
    );
    const histogram = useMemo(
        () => buildHistogram(rasterData ?? [], domain, HISTOGRAM_BINS, resolvedXScale),
        [domain, rasterData, resolvedXScale],
    );
    const resolvedYScale = useMemo(
        () => (yScaleMode === "auto" ? resolveAutoYScale(histogram.counts) : yScaleMode),
        [histogram.counts, yScaleMode],
    );

    const layers = useMemo<HistogramLayer[]>(() => {
        const maxCount = Math.max(1, ...histogram.counts);
        const highlightData = new Array(HISTOGRAM_BINS).fill(0);
        if (
            Number.isFinite(pixelValue) &&
            domain[0] !== domain[1] &&
            pixelValue >= domain[0] &&
            pixelValue <= domain[1]
        ) {
            const highlightIndex = Math.max(
                0,
                histogram.edges.findIndex((edge, i) =>
                    i < histogram.edges.length - 1 &&
                    pixelValue >= edge &&
                    pixelValue <= histogram.edges[i + 1],
                ),
            );
            highlightData[highlightIndex] = maxCount;
        }
        return [
            {
                id: `channel-${index}`,
                data: histogram.counts,
                color: `rgba(${color[0]}, ${color[1]}, ${color[2]}, 0.55)`,
                variant: "line",
            },
            {
                id: `channel-${index}-pixel`,
                data: highlightData,
                color: "rgba(255, 255, 255, 0.95)",
                variant: "markers",
                widthFactor: 0.35,
                hidden: highlightData.every((value) => value === 0),
            },
        ];
    }, [color, domain, histogram.counts, histogram.edges, index, pixelValue]);

    const isHistogramLoading =
        !currentDomain || !contrastLimits[index] || isChannelLoading[index];

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
                        onToggleX: () =>
                            setXScaleMode((mode) => toggleScaleMode(mode, resolvedXScale)),
                        onToggleY: () =>
                            setYScaleMode((mode) => toggleScaleMode(mode, resolvedYScale)),
                    }}
                    markers={[
                        {
                            id: `pixel-value-${index}`,
                            value: pixelValue,
                            color: "rgba(255, 255, 255, 0.95)",
                            hidden:
                                !Number.isFinite(pixelValue) ||
                                pixelValue < domain[0] ||
                                pixelValue > domain[1],
                        },
                    ]}
                />
            )}
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
        <Accordion
            disableGutters
            elevation={0}
            defaultExpanded={false}
            sx={{
                gridColumn: "1 / -1",
                border: "1px solid hsl(var(--border))",
                borderRadius: "6px",
                backgroundColor: "hsl(var(--muted) / 0.35)",
                "&:before": { display: "none" },
            }}
        >
            <AccordionSummary
                expandIcon={<ExpandMoreIcon fontSize="small" />}
                sx={{
                    minHeight: "34px",
                    px: 1,
                    "& .MuiAccordionSummary-content": {
                        my: 0.5,
                    },
                }}
            >
                <span className="text-[10px] font-medium uppercase tracking-[0.14em] text-[hsl(var(--muted-foreground))]">
                    Tone controls
                </span>
            </AccordionSummary>
            <AccordionDetails sx={{ px: 1, pb: 1, pt: 0 }}>
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
                                contrast[index] = v as number;
                                const newContrast = [...contrast];
                                channelsStore.setState({ contrast: newContrast });
                            }}
                            onClick={(e) => {
                                if (e.detail === 2) {
                                    contrast[index] = 0.5;
                                    const newContrast = [...contrast];
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
                                brightness[index] = v as number;
                                const newbrightness = [...brightness];
                                channelsStore.setState({ brightness: newbrightness });
                            }}
                            onClick={(e) => {
                                if (e.detail === 2) {
                                    brightness[index] = 0.5;
                                    const newbrightness = [...brightness];
                                    channelsStore.setState({ brightness: newbrightness });
                                }
                            }}
                        />
                    </div>
                </div>
            </AccordionDetails>
        </Accordion>
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

    return (
        <>
            <div
                className="relative grid grid-cols-[minmax(10rem,0.95fr)_minmax(0,1.15fr)] gap-2 rounded-md border border-[hsl(var(--border))] bg-[hsl(var(--background))] p-2 shadow-sm"
                onMouseEnter={() => setIsHovered(true)}
                onMouseLeave={() => setIsHovered(false)}
            >
                <div
                    className="min-h-[94px] rounded-md border border-[hsl(var(--border))] bg-[hsl(var(--background))] p-2"
                    style={{
                        boxShadow: `inset 3px 0 0 rgb(${color[0]}, ${color[1]}, ${color[2]})`,
                    }}
                >
                    <div className="mb-1 flex items-center justify-between gap-2 text-[10px] font-medium uppercase tracking-[0.14em] text-[hsl(var(--muted-foreground))]">
                        <span>Channel</span>
                        <div className="flex items-center gap-1.5">
                            <Checkbox
                                checked={channelVisible}
                                disabled={isChannelLoading[index]}
                                sx={{
                                    padding: "2px",
                                    color: "hsl(var(--muted-foreground))",
                                    "&.Mui-checked": {
                                        color: `rgb(${color[0]}, ${color[1]}, ${color[2]})`,
                                    },
                                }}
                                onChange={() => {
                                    const channelsVisible =
                                        channelsStore.getState().channelsVisible;
                                    channelsVisible[index] = !channelVisible;
                                    const visible = [...channelsVisible];
                                    channelsStore.setState({ channelsVisible: visible });
                                }}
                            />
                            <div className="h-7 min-w-8 overflow-hidden rounded border border-[hsl(var(--border))] bg-[hsl(var(--background))] shadow-sm">
                                <PopoverPicker
                                    color={color}
                                    onChange={(c) => {
                                        const colors = channelsStore.getState().colors;
                                        colors[index] = c;
                                        const newColors = [...colors];
                                        channelsStore.setState({ colors: newColors });
                                    }}
                                />
                            </div>
                        </div>
                    </div>
                    <ChannelChooser index={index} />
                </div>
                <div className="min-w-0 self-stretch">
                    <ChannelHistogram index={index} />
                </div>
                <IconButton
                    onClick={() => {
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
                        visibility: isHovered ? "visible" : "hidden",
                        transition: "opacity 0.3s, visibility 0.3s",
                        "&:focus": {
                            visibility: "visible",
                            opacity: 1,
                        },
                    }}
                >
                    <HighlightOffIcon fontSize="small" />
                </IconButton>
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
            className="rounded-md border border-[hsl(var(--border))] bg-[hsl(var(--background))] px-3 py-1.5 text-sm font-medium text-[hsl(var(--foreground))] shadow-sm transition hover:bg-[hsl(var(--muted)/0.45)] disabled:cursor-not-allowed disabled:opacity-40"
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
        <div className="w-full space-y-2 bg-[hsl(var(--background))] p-2">
            {ids.map((id, i) => (
                <ChannelController key={id} index={i} />
            ))}
            <AddChannel />
        </div>
    );
};
