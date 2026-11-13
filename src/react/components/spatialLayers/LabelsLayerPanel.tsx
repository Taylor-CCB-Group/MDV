import { Checkbox, FormControlLabel, Slider, Typography } from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";

import type DataStore from "@/datastore/DataStore";
import { mdvFieldSpecsOf, type WithMdvFieldSpecs } from "@/react/spatialdata/field_spec_projection";
import type { TableAssociation } from "@/react/spatialdata/table_association";
import { ColorFields, type Rgb } from "./ColorFields";
import TableAssociationControls, {
    FillColorByColumnControl,
    type LayerFillColorByColumn,
} from "./TableAssociationControls";

type LabelsLayerConfig = WithMdvFieldSpecs<
    Extract<LayerConfig, { type: "labels" }> & {
        fillColorByColumn?: LayerFillColorByColumn;
    }
>;
type LabelsChannels = NonNullable<LabelsLayerConfig["channels"]>;

type Props = {
    config: LabelsLayerConfig;
    updateLayer: (updates: Partial<LabelsLayerConfig>) => void;
    association: TableAssociation;
    dataStore: DataStore;
};

/**
 * The labels layer's own defaults (`LabelsLayer` in `@spatialdata/layers`), for showing
 * a value nobody has set. The element's metadata can override them; it rarely does.
 */
const DEFAULT_LABEL_COLOR: Rgb = [255, 255, 255];
const LABEL_DEFAULTS = {
    filled: true,
    fillOpacity: 0.18,
    outlineOpacity: 0.95,
    strokeWidth: 1.5,
};

/**
 * `values` with its first entry replaced. The layer reads each channel array whole when it
 * is non-empty, so the other channels' entries are kept rather than dropped.
 */
function withFirst<T>(values: T[] | undefined, value: T): T[] {
    const next = values ? [...values] : [];
    next[0] = value;
    return next;
}

function LabeledSlider({
    label,
    value,
    min,
    max,
    step,
    onChange,
}: {
    label: string;
    value: number;
    min: number;
    max: number;
    step: number;
    onChange: (value: number) => void;
}) {
    return (
        <div className="flex items-center gap-3">
            <Typography variant="caption" className="w-28 flex-shrink-0">
                {label}
            </Typography>
            <Slider
                size="small"
                min={min}
                max={max}
                step={step}
                value={value}
                aria-label={label}
                onChange={(_, next) => {
                    if (typeof next === "number") onChange(next);
                }}
            />
        </div>
    );
}

export default function LabelsLayerPanel({ config, updateLayer, association, dataStore }: Props) {
    const specs = mdvFieldSpecsOf(config);
    const channels = config.channels ?? {};
    // Only the fields this patch changes: `patchRenderStackEntryPropsInPlace` merges
    // channel fields into the saved ones rather than replacing the record.
    const patchChannels = (patch: Partial<LabelsChannels>) => updateLayer({ channels: patch });
    const color = channels.colors?.[0] ?? DEFAULT_LABEL_COLOR;
    const filled = channels.channelsFilled?.[0] ?? LABEL_DEFAULTS.filled;

    return (
        <div className="grid gap-3">
            <TableAssociationControls
                association={association}
                dataStore={dataStore}
                tooltipFields={config.tooltipFields ?? []}
                tooltipFieldsSpec={specs?.tooltipFields}
                onTooltipFieldsChange={(next, spec) =>
                    updateLayer({
                        tooltipFields: next,
                        mdvFieldSpecs: { ...specs, tooltipFields: spec },
                    })
                }
            />
            <FillColorByColumnControl
                association={association}
                dataStore={dataStore}
                fillColorByColumn={config.fillColorByColumn}
                fillColorByColumnSpec={specs?.fillColorByColumn}
                onChange={(next, spec) =>
                    updateLayer({
                        fillColorByColumn: next,
                        mdvFieldSpecs: { ...specs, fillColorByColumn: spec },
                    })
                }
            />
            <ColorFields
                label="Colour"
                rgb={color}
                caption={config.fillColorByColumn ? "Coloured by column; outlines and opacity still apply." : undefined}
                onChange={(rgb) => patchChannels({ colors: withFirst(channels.colors, rgb) })}
            />
            <FormControlLabel
                control={
                    <Checkbox
                        size="small"
                        checked={filled}
                        onChange={(e) =>
                            patchChannels({ channelsFilled: withFirst(channels.channelsFilled, e.target.checked) })
                        }
                    />
                }
                label={<Typography variant="caption">Filled</Typography>}
            />
            {filled && (
                <LabeledSlider
                    label="Fill opacity"
                    min={0}
                    max={1}
                    step={0.01}
                    value={channels.channelOpacities?.[0] ?? LABEL_DEFAULTS.fillOpacity}
                    onChange={(value) =>
                        patchChannels({ channelOpacities: withFirst(channels.channelOpacities, value) })
                    }
                />
            )}
            <LabeledSlider
                label="Outline opacity"
                min={0}
                max={1}
                step={0.01}
                value={channels.channelOutlineOpacities?.[0] ?? LABEL_DEFAULTS.outlineOpacity}
                onChange={(value) =>
                    patchChannels({ channelOutlineOpacities: withFirst(channels.channelOutlineOpacities, value) })
                }
            />
            <LabeledSlider
                label="Outline width"
                min={0}
                max={5}
                step={0.1}
                value={channels.channelStrokeWidths?.[0] ?? LABEL_DEFAULTS.strokeWidth}
                onChange={(value) =>
                    patchChannels({ channelStrokeWidths: withFirst(channels.channelStrokeWidths, value) })
                }
            />
        </div>
    );
}
