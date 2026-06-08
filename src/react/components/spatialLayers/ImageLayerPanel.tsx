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

type ImageLayerConfig = Extract<LayerConfig, { type: "image" }>;
import { useId, useMemo } from "react";
import { PopoverPicker } from "../ColorPicker";
import type { ImageLayerSource } from "@/react/spatial_layer_stack";

type ChannelConfig = NonNullable<ImageLayerConfig["channels"]>;

type Props = {
    layerId: string;
    config: ImageLayerConfig;
    imageSource?: ImageLayerSource;
    loaderDefaults?: {
        colors?: [number, number, number][];
        contrastLimits?: [number, number][];
        channelsVisible?: boolean[];
        selections?: ChannelConfig["selections"];
    };
    channelNames?: string[];
    updateLayer: (updates: Partial<LayerConfig>) => void;
};

function withChannelValue<T>(values: T[] | undefined, index: number, value: T, fallback: T): T[] {
    const next = [...(values ?? [])];
    while (next.length <= index) {
        next.push(fallback);
    }
    next[index] = value;
    return next;
}

export default function ImageLayerPanel({
    config,
    imageSource = "spatial",
    loaderDefaults,
    channelNames = [],
    updateLayer,
}: Props) {
    const channels = config.channels ?? {};
    const colors = channels.colors ?? loaderDefaults?.colors ?? [[255, 0, 0]];
    const contrastLimits =
        channels.contrastLimits ?? loaderDefaults?.contrastLimits ?? [[0, 255]];
    const channelsVisible =
        channels.channelsVisible ?? loaderDefaults?.channelsVisible ?? colors.map(() => true);
    const selections = channels.selections ?? loaderDefaults?.selections ?? [{ c: 0 }];

    const channelCount = Math.max(colors.length, contrastLimits.length, channelsVisible.length, 1);
    const channelIndexes = useMemo(
        () => Array.from({ length: channelCount }, (_, index) => index),
        [channelCount],
    );

    const updateChannels = (patch: Partial<ChannelConfig>) => {
        updateLayer({
            channels: {
                ...channels,
                ...patch,
            },
        });
    };

    return (
        <div className="space-y-3">
            <Typography variant="caption" color="text.secondary">
                Image source: {imageSource === "spatial" ? "SpatialData zarr" : "OME-TIFF"}
            </Typography>
            {channelIndexes.map((index) => (
                <ChannelRow
                    key={`${config.id}-channel-${index}`}
                    index={index}
                    name={channelNames[selections[index]?.c ?? index] ?? `Channel ${index + 1}`}
                    color={colors[index] ?? [255, 0, 0]}
                    contrast={contrastLimits[index] ?? [0, 255]}
                    visible={channelsVisible[index] ?? true}
                    channelOptions={channelNames}
                    selectedChannel={selections[index]?.c ?? 0}
                    onVisibleChange={(visible) => {
                        updateChannels({
                            channelsVisible: withChannelValue(channelsVisible, index, visible, true),
                        });
                    }}
                    onColorChange={(color) => {
                        updateChannels({
                            colors: withChannelValue(colors, index, color, [255, 0, 0]),
                        });
                    }}
                    onContrastChange={(limits) => {
                        updateChannels({
                            contrastLimits: withChannelValue(contrastLimits, index, limits, [0, 255]),
                        });
                    }}
                    onChannelSelect={(channelIndex) => {
                        const nextSelections = [...selections];
                        while (nextSelections.length <= index) {
                            nextSelections.push({ c: 0 });
                        }
                        nextSelections[index] = { ...nextSelections[index], c: channelIndex };
                        updateChannels({ selections: nextSelections });
                    }}
                />
            ))}
        </div>
    );
}

function ChannelRow({
    index,
    name,
    color,
    contrast,
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
    contrast: [number, number];
    visible: boolean;
    channelOptions: string[];
    selectedChannel: number;
    onVisibleChange: (visible: boolean) => void;
    onColorChange: (color: [number, number, number]) => void;
    onContrastChange: (limits: [number, number]) => void;
    onChannelSelect: (channelIndex: number) => void;
}) {
    const labelId = useId();
    const [low, high] = contrast;

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
                        min={0}
                        max={high}
                        onChange={(_, value) => {
                            if (typeof value === "number") {
                                onContrastChange([value, high]);
                            }
                        }}
                    />
                </div>
                <div className="flex items-center gap-2">
                    <span className="w-10 text-xs">High</span>
                    <Slider
                        size="small"
                        value={high}
                        min={low}
                        max={Math.max(high * 2, low + 1)}
                        onChange={(_, value) => {
                            if (typeof value === "number") {
                                onContrastChange([low, value]);
                            }
                        }}
                    />
                </div>
            </div>
        </div>
    );
}
