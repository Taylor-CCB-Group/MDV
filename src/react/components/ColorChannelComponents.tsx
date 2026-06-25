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
    useViewerStoreApi,
} from "./avivatorish/state";
import { getSingleSelectionStats } from "./avivatorish/utils";
import { COLOR_PALLETE } from "./avivatorish/constants";
import { MAX_CHANNELS } from "@vivjs/constants";
import { useSpatialImagePanelContext } from "./spatialLayers/ImageLayerPanel";

const DEFAULT_BRIGHTNESS_CONTRAST = 0.5;

function createStableChannelId(layerId: string): string {
    const suffix =
        typeof crypto !== "undefined" && typeof crypto.randomUUID === "function"
            ? crypto.randomUUID().slice(0, 8)
            : Math.random().toString(36).slice(2, 10);
    return `${layerId}-ch-${suffix}`;
}

function withChannelValue(values: number[], index: number, value: number) {
    return Array.from({ length: Math.max(values.length, index + 1) }, (_, i) =>
        i === index ? value : values[i] ?? DEFAULT_BRIGHTNESS_CONTRAST,
    );
}

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
    const spatial = useSpatialImagePanelContext();
    const metadata = useMetadata();
    const { selections, setPropertiesForChannel } = useChannelsStore(
        ({ selections, setPropertiesForChannel }) => ({
            selections,
            setPropertiesForChannel,
        }),
        shallow,
    );
    const loader = useLoader();
    const { setIsChannelLoading, isChannelLoading, removeIsChannelLoading, use3d, channelOptions } = useViewerStore(
        ({ setIsChannelLoading, isChannelLoading, removeIsChannelLoading, use3d, channelOptions }) => ({
            setIsChannelLoading,
            isChannelLoading,
            removeIsChannelLoading,
            use3d,
            channelOptions,
        }),
        shallow,
    );
    const channels =
        spatial?.channelNames ??
        metadata?.Pixels.Channels.map((c, i) => c.Name ?? `Channel ${i + 1}`) ??
        channelOptions ??
        selections.map((_, i) => `Channel ${i + 1}`);

    if (spatial) {
        return (
            <select
                disabled={isChannelLoading[index]}
                value={spatial.selections[index]?.c ?? 0}
                className="w-full rounded-md border border-[hsl(var(--border))] bg-[hsl(var(--background))] px-2 py-1.5 text-sm text-[hsl(var(--foreground))] shadow-sm outline-none transition focus:border-[hsl(var(--ring))] disabled:cursor-wait disabled:opacity-60"
                onChange={(e) => {
                    const c = Number.parseInt(e.target.value, 10);
                    if (Number.isNaN(c)) return;
                    const nextSelections = spatial.selections.map((selection, i) =>
                        i === index ? { ...selection, c } : selection,
                    );
                    spatial.setChannels({ selections: nextSelections });
                }}
            >
                {channels.map((c, i) => (
                    <option key={`${i}_${c}`} value={i}>
                        {c}
                    </option>
                ))}
            </select>
        );
    }

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
                        const data = Array.isArray(loader) ? loader[loader.length - 1] : loader;
                        if (!data?.getRaster) {
                            setPropertiesForChannel(index, { selections: selection });
                            return;
                        }
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
    const spatial = useSpatialImagePanelContext();
    const contrastLimits = useChannelsStore((state) => state.contrastLimits);
    const storeColor = useChannelsStore((state) => state.colors[index] ?? [37, 99, 235]);
    // Spatial mode: color is canonical channel config, read from the panel context.
    const color = spatial ? (spatial.colors[index] ?? storeColor) : storeColor;
    const currentDomain = useChannelsStore((state) => state.domains[index]);
    const domain = currentDomain ?? ([0, 1] as Range);
    const rasterData = useChannelsStore((state) => state.raster[index]?.data);
    const pixelValue = useViewerStore((state) => state.pixelValues[index]);
    const isChannelLoading = useViewerStore((state) => state.isChannelLoading);
    const channelsStore = useChannelsStoreApi();
    const spatialRef = useRef(spatial);
    spatialRef.current = spatial;
    const dark = useTheme() === "dark";
    const limits = spatial
        ? (spatial.contrastLimits[index] ?? domain)
        : (contrastLimits[index] ?? domain);
    const domainRef = useRef(domain);
    const limitsRef = useRef(limits);
    const [xScaleMode, setXScaleMode] = useState<ScaleMode>("auto");
    const [yScaleMode, setYScaleMode] = useState<ScaleMode>("auto");

    useEffect(() => {
        domainRef.current = domain;
        limitsRef.current = limits;
    }, [domain, limits]);

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

    const isHistogramLoading =
        !currentDomain ||
        !(spatial ? spatial.contrastLimits[index] : contrastLimits[index]) ||
        isChannelLoading[index];

    const handleBrushValue = useCallback(
        (value: Range | null) => {
            const currentDomain = domainRef.current;
            const currentLimits = limitsRef.current;
            const nextValue = value ? sortRange(clampRange(value, currentDomain)) : currentDomain;
            if (nextValue.some((item) => Number.isNaN(item)) || rangesEqual(nextValue, currentLimits)) return;
            limitsRef.current = nextValue;
            const panel = spatialRef.current;
            if (panel) {
                const nextContrastLimits = panel.contrastLimits.map((limits, i) =>
                    i === index ? nextValue : limits,
                );
                panel.setChannels({ contrastLimits: nextContrastLimits });
                return;
            }
            const nextContrastLimits = [...channelsStore.getState().contrastLimits];
            nextContrastLimits[index] = nextValue;
            channelsStore.setState({ contrastLimits: nextContrastLimits });
        },
        [channelsStore, index],
    );

    const brush = useMemo(
        () => ({
            value: rangesEqual(limits, domain) ? null : limits,
            setValue: handleBrushValue,
            minMax: domain,
        }),
        [domain, handleBrushValue, limits],
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
    const spatial = useSpatialImagePanelContext();
    const storeTone = useChannelsStore(({ contrast, brightness }) => ({ contrast, brightness }));
    // Spatial mode: tone is canonical (entry `vivLayerProps`), read from the panel context.
    const contrast = spatial ? spatial.contrast : storeTone.contrast;
    const brightness = spatial ? spatial.brightness : storeTone.brightness;
    const channelsStore = useChannelsStoreApi();
    const isChannelLoading = useViewerStore((state) => state.isChannelLoading);
    const contrastValue = contrast[index] ?? DEFAULT_BRIGHTNESS_CONTRAST;
    const brightnessValue = brightness[index] ?? DEFAULT_BRIGHTNESS_CONTRAST;
    return (
        <div className="grid gap-x-4 gap-y-2 sm:grid-cols-2">
            <div className="flex items-center gap-2">
                <span className="w-16 shrink-0 text-[10px] font-medium uppercase tracking-[0.14em] text-[hsl(var(--muted-foreground))]">
                    Contrast
                </span>
                <Slider
                    size="small"
                    disabled={isChannelLoading[index]}
                    value={contrastValue}
                    min={0.1}
                    max={0.9}
                    step={0.01}
                    onChange={(_, v) => {
                        if (isArray(v)) return;
                        if (spatial) {
                            spatial.patchToneAtIndex(index, "contrast", v);
                            return;
                        }
                        channelsStore.setState({ contrast: withChannelValue(contrast, index, v) });
                    }}
                    onClick={(e) => {
                        if (e.detail === 2) {
                            if (spatial) {
                                spatial.patchToneAtIndex(index, "contrast", DEFAULT_BRIGHTNESS_CONTRAST);
                                return;
                            }
                            channelsStore.setState({
                                contrast: withChannelValue(contrast, index, DEFAULT_BRIGHTNESS_CONTRAST),
                            });
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
                    value={brightnessValue}
                    min={0.01}
                    max={0.99}
                    step={0.01}
                    onChange={(_, v) => {
                        if (isArray(v)) return;
                        if (spatial) {
                            spatial.patchToneAtIndex(index, "brightness", v);
                            return;
                        }
                        channelsStore.setState({ brightness: withChannelValue(brightness, index, v) });
                    }}
                    onClick={(e) => {
                        if (e.detail === 2) {
                            if (spatial) {
                                spatial.patchToneAtIndex(index, "brightness", DEFAULT_BRIGHTNESS_CONTRAST);
                                return;
                            }
                            channelsStore.setState({
                                brightness: withChannelValue(brightness, index, DEFAULT_BRIGHTNESS_CONTRAST),
                            });
                        }
                    }}
                />
            </div>
        </div>
    );
};

const ChannelController = ({ index }: { index: number }) => {
    const spatial = useSpatialImagePanelContext();
    const storeColor = useChannelsStore((state) => state.colors[index]);
    const storeChannelVisible = useChannelsStore((state) => state.channelsVisible[index]);
    // Spatial mode: color + visibility are canonical channel config, read from the panel context.
    const color = spatial ? spatial.colors[index] : storeColor;
    const channelVisible = spatial ? spatial.channelsVisible[index] : storeChannelVisible;
    const channelIds = useChannelsStore((state) => state.ids);
    const removeChannel = useChannelsStore((state) => state.removeChannel);
    const isChannelLoading = useViewerStore((state) => state.isChannelLoading);
    const channelsStore = useChannelsStoreApi();
    const [isHovered, setIsHovered] = useState(false);

    if (!color) return null;

    const hasPendingLoads = isChannelLoading.some(Boolean);
    const canRemoveChannel = spatial ? spatial.channelIds.length > 1 : channelIds.length > 1;

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
                                        if (spatial) {
                                            const visible = spatial.channelsVisible.map((value, i) =>
                                                i === index ? !channelVisible : value,
                                            );
                                            spatial.setChannels({ channelsVisible: visible });
                                            return;
                                        }
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
                                            if (spatial) {
                                                const colors = spatial.colors.map((value, i) =>
                                                    i === index ? c : value,
                                                );
                                                spatial.setChannels({ colors });
                                                return;
                                            }
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
                    disabled={hasPendingLoads || !canRemoveChannel}
                    onClick={(event) => {
                        event.stopPropagation();
                        if (spatial) {
                            spatial.removeChannel(index);
                            return;
                        }
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
    const spatial = useSpatialImagePanelContext();
    const loader = useLoader();
    const isAddingChannelRef = useRef(false);
    const [isAddingChannel, setIsAddingChannel] = useState(false);
    // const { labels } = loader[0];
    // const channelsStore = useChannelsStoreApi();
    const { selections, domains: channelDomains, contrastLimits: channelContrastLimits, setPropertiesForChannel } = useChannelsStore(
        ({ selections, domains, contrastLimits, setPropertiesForChannel }) => ({
            selections,
            domains,
            contrastLimits,
            setPropertiesForChannel,
        }),
        shallow,
    );
    const canAddChannel = spatial
        ? spatial.channelIds.length < MAX_CHANNELS
        : selections.length < MAX_CHANNELS;
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
    const viewerStore = useViewerStoreApi();
    const canLoadStats = Boolean((Array.isArray(loader) ? loader[loader.length - 1] : loader)?.getRaster);
    return (
        <button
            type="button"
            className="rounded-md border border-[hsl(var(--border))] bg-[hsl(var(--background))] px-3 py-1.5 text-sm font-medium text-[hsl(var(--foreground))] shadow-sm transition hover:bg-[hsl(var(--muted)/0.45)] disabled:cursor-not-allowed disabled:opacity-40"
            disabled={!canAddChannel || isAddingChannel}
            onClick={async () => {
                if (spatial) {
                    if (isAddingChannelRef.current) return;
                    isAddingChannelRef.current = true;
                    setIsAddingChannel(true);
                    try {
                        spatial.addChannelWithSelection();
                    } finally {
                        isAddingChannelRef.current = false;
                        setIsAddingChannel(false);
                    }
                    return;
                }
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
                        c: index,
                    };
                    const channelId = createStableChannelId("channel");
                    addChannel({
                        selections: selection,
                        ids: channelId,
                        channelsVisible: true,
                        colors: COLOR_PALLETE[index % COLOR_PALLETE.length],
                        domains: channelDomains[index] ?? channelDomains[0] ?? [0, 255],
                        contrastLimits: channelContrastLimits[index] ?? channelContrastLimits[0] ?? [0, 255],
                    });
                    viewerStore.setState((state) => ({
                        channelOptions: Array.from({ length: index + 1 }, (_, i) =>
                            state.channelOptions[i] ?? `Channel ${i + 1}`,
                        ),
                        isChannelLoading: Array.from({ length: index + 1 }, (_, i) =>
                            i === index ? true : (state.isChannelLoading[i] ?? false),
                        ),
                        pixelValues: Array.from({ length: index + 1 }, (_, i) =>
                            state.pixelValues[i] ?? Number.NaN,
                        ),
                    }));
                    didAddChannel = true;
                    if (!canLoadStats) {
                        setIsChannelLoading(index, false);
                        return;
                    }
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

export const VivChannelList = () => {
    const spatial = useSpatialImagePanelContext();
    const storeIds = useChannelsStore(({ ids }) => ids);
    // Spatial mode: the channel list is driven by canonical channel ids.
    const ids = spatial ? spatial.channelIds : storeIds;
    return (
        <div className="w-full space-y-2 bg-[hsl(var(--background))] p-2">
            {ids.map((id, i) => (
                <ChannelController key={id} index={i} />
            ))}
            <AddChannel />
        </div>
    );
};

export const Test = VivChannelList;
