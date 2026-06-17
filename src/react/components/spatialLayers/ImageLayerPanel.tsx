import {
    Checkbox,
    FormControl,
    InputLabel,
    MenuItem,
    Select,
    Slider,
    Typography,
} from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";
import { useEffect, useId, useMemo, useState } from "react";
import { useDebounce } from "use-debounce";
import { PopoverPicker } from "../ColorPicker";
import type { ImageLayerSource } from "@/react/spatial_layer_stack";

type ImageLayerConfig = Extract<LayerConfig, { type: "image" }>;
type ChannelConfig = NonNullable<ImageLayerConfig["channels"]>;
type Range = [number, number];

type LoaderDefaults = {
    colors?: [number, number, number][];
    contrastLimits?: [number, number][];
    channelsVisible?: boolean[];
    selections?: ChannelConfig["selections"];
};

type Props = {
    layerId: string;
    config: ImageLayerConfig;
    imageSource?: ImageLayerSource;
    loaderDefaults?: LoaderDefaults;
    channelNames?: string[];
    updateLayer: (updates: Partial<LayerConfig>) => void;
    patchLayer: (updates: Partial<LayerConfig>) => void;
};

function copyRange(value: [number, number] | undefined, fallback: Range): Range {
    if (!value) return fallback;
    return [value[0], value[1]];
}

function withChannelValue<T>(values: T[] | undefined, index: number, value: T, fallback: T): T[] {
    const next = [...(values ?? [])];
    while (next.length <= index) {
        next.push(fallback);
    }
    next[index] = value;
    return next;
}

function clampContrast(value: Range, domain: Range): Range {
    const [domainMin, domainMax] = domain;
    let [low, high] = value;
    low = Math.max(domainMin, Math.min(low, domainMax));
    high = Math.max(domainMin, Math.min(high, domainMax));
    if (low >= high) {
        if (high < domainMax) {
            low = Math.max(domainMin, high - 1);
        } else {
            high = Math.min(domainMax, low + 1);
        }
    }
    return [low, high];
}

function getChannelContrast(
    channels: ChannelConfig,
    loaderDefaults: LoaderDefaults | undefined,
    index: number,
    domain: Range,
): Range {
    const fromConfig = channels.contrastLimits?.[index];
    if (fromConfig) return copyRange(fromConfig, domain);
    const fromLoader = loaderDefaults?.contrastLimits?.[index];
    if (fromLoader) return copyRange(fromLoader, domain);
    return domain;
}

export default function ImageLayerPanel({
    config,
    imageSource = "spatial",
    loaderDefaults,
    channelNames = [],
    updateLayer,
    patchLayer,
}: Props) {
    const channels = config.channels ?? {};
    const colors = channels.colors ?? loaderDefaults?.colors ?? [[255, 0, 0]];
    const channelsVisible =
        channels.channelsVisible ?? loaderDefaults?.channelsVisible ?? colors.map(() => true);
    const selections = channels.selections ?? loaderDefaults?.selections ?? [{ c: 0 }];

    const channelCount = Math.max(
        colors.length,
        loaderDefaults?.contrastLimits?.length ?? 0,
        channels.contrastLimits?.length ?? 0,
        channelsVisible.length,
        1,
    );
    const channelIndexes = useMemo(
        () => Array.from({ length: channelCount }, (_, index) => index),
        [channelCount],
    );

    const buildChannels = (patch: Partial<ChannelConfig>): ChannelConfig => ({
        colors,
        contrastLimits: channels.contrastLimits,
        channelsVisible,
        selections,
        ...channels,
        ...patch,
    });

    const applyChannelPatch = (patch: Partial<ChannelConfig>, persist = false) => {
        const nextChannels = buildChannels(patch);
        patchLayer({ channels: nextChannels });
        if (persist) {
            updateLayer({ channels: nextChannels });
        }
    };

    return (
        <div className="space-y-3">
            <Typography variant="caption" color="text.secondary">
                Image source: {imageSource === "spatial" ? "SpatialData zarr" : "OME-TIFF"}
            </Typography>
            {channelIndexes.map((index) => {
                const domain = copyRange(
                    loaderDefaults?.contrastLimits?.[index],
                    [0, 255],
                );
                return (
                    <ChannelRow
                        key={`${config.id}-channel-${index}`}
                        index={index}
                        name={channelNames[selections[index]?.c ?? index] ?? `Channel ${index + 1}`}
                        color={colors[index] ?? [255, 0, 0]}
                        contrast={getChannelContrast(channels, loaderDefaults, index, domain)}
                        domain={domain}
                        visible={channelsVisible[index] ?? true}
                        channelOptions={channelNames}
                        selectedChannel={selections[index]?.c ?? 0}
                        onVisibleChange={(visible) => {
                            applyChannelPatch({
                                channelsVisible: withChannelValue(channelsVisible, index, visible, true),
                            }, true);
                        }}
                        onColorChange={(color) => {
                            applyChannelPatch({
                                colors: withChannelValue(colors, index, color, [255, 0, 0]),
                            }, true);
                        }}
                        onContrastChange={(limits, persist) => {
                            applyChannelPatch({
                                contrastLimits: withChannelValue(
                                    channels.contrastLimits ?? loaderDefaults?.contrastLimits,
                                    index,
                                    limits,
                                    domain,
                                ),
                            }, persist);
                        }}
                        onChannelSelect={(channelIndex) => {
                            const nextSelections = [...selections];
                            while (nextSelections.length <= index) {
                                nextSelections.push({ c: 0 });
                            }
                            nextSelections[index] = { ...nextSelections[index], c: channelIndex };
                            applyChannelPatch({ selections: nextSelections }, true);
                        }}
                    />
                );
            })}
        </div>
    );
}

function ChannelRow({
    index,
    name,
    color,
    contrast,
    domain,
    visible,
    channelOptions,
    selectedChannel,
    onVisibleChange,
    onColorChange,
    onContrastChange,
    onChannelSelect,
}: {
    index: number;
    name: string;
    color: [number, number, number];
    contrast: Range;
    domain: Range;
    visible: boolean;
    channelOptions: string[];
    selectedChannel: number;
    onVisibleChange: (visible: boolean) => void;
    onColorChange: (color: [number, number, number]) => void;
    onContrastChange: (limits: Range, persist: boolean) => void;
    onChannelSelect: (channelIndex: number) => void;
}) {
    const labelId = useId();
    const [domainMin, domainMax] = domain;
    const [liveContrast, setLiveContrast] = useState<Range | null>(null);
    const [debouncedContrast] = useDebounce(liveContrast, 50);

    useEffect(() => {
        setLiveContrast(null);
    }, [contrast[0], contrast[1]]);

    useEffect(() => {
        if (!debouncedContrast) return;
        if (debouncedContrast[0] === contrast[0] && debouncedContrast[1] === contrast[1]) {
            return;
        }
        onContrastChange(debouncedContrast, true);
    }, [contrast, debouncedContrast, onContrastChange]);

    const [low, high] = liveContrast ?? contrast;
    const sliderStep = Math.max(1, Math.round((domainMax - domainMin) / 500));

    return (
        <div className="rounded-md border border-[hsl(var(--border))] p-2">
            <div className="mb-2 flex items-center justify-between gap-2">
                <Typography variant="subtitle2">
                    {name}
                </Typography>
                <Checkbox
                    size="small"
                    checked={visible}
                    onChange={(event) => onVisibleChange(event.target.checked)}
                />
            </div>
            {channelOptions.length > 0 && (
                <FormControl fullWidth size="small" variant="standard" className="mb-2">
                    <InputLabel id={labelId}>Channel</InputLabel>
                    <Select
                        labelId={labelId}
                        value={selectedChannel}
                        onChange={(event) => onChannelSelect(Number(event.target.value))}
                    >
                        {channelOptions.map((option, optionIndex) => (
                            <MenuItem key={`${index}-${option}`} value={optionIndex}>
                                {option}
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>
            )}
            <div className="mb-2 h-8 w-12 overflow-hidden rounded border border-[hsl(var(--border))]">
                <PopoverPicker color={color} onChange={onColorChange} />
            </div>
            <div className="grid gap-2">
                <div className="flex items-center gap-2">
                    <span className="w-10 text-xs">Low</span>
                    <Slider
                        size="small"
                        value={low}
                        min={domainMin}
                        max={domainMax}
                        step={sliderStep}
                        onChange={(_, value) => {
                            if (typeof value !== "number") return;
                            const next = clampContrast([value, high], domain);
                            setLiveContrast(next);
                            onContrastChange(next, false);
                        }}
                    />
                </div>
                <div className="flex items-center gap-2">
                    <span className="w-10 text-xs">High</span>
                    <Slider
                        size="small"
                        value={high}
                        min={domainMin}
                        max={domainMax}
                        step={sliderStep}
                        onChange={(_, value) => {
                            if (typeof value !== "number") return;
                            const next = clampContrast([low, value], domain);
                            setLiveContrast(next);
                            onContrastChange(next, false);
                        }}
                    />
                </div>
            </div>
        </div>
    );
}
