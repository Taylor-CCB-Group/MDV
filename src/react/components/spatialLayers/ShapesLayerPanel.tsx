import {
    Autocomplete,
    Chip,
    MenuItem,
    Select,
    Slider,
    TextField,
    Typography,
} from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";

type ShapesLayerConfig = Extract<LayerConfig, { type: "shapes" }>;
import { useMemo } from "react";
import type { TableAssociation } from "@/react/spatialdata/table_association";

type Props = {
    config: ShapesLayerConfig;
    updateLayer: (updates: Partial<ShapesLayerConfig>) => void;
    association: TableAssociation;
    availableFields: string[];
    /** Chart-level color_by column, for fill-by-column option list.
     * nb this needs review, not working and also general background issues with colorBy type.
     */
    chartColorBy?: string;
};

function ColorFields({
    label,
    value,
    onChange,
}: {
    label: string;
    value: [number, number, number, number];
    onChange: (next: [number, number, number, number]) => void;
}) {
    return (
        <div className="grid gap-2">
            <Typography variant="caption">{label}</Typography>
            <div className="grid grid-cols-4 gap-2">
                {(["R", "G", "B", "A"] as const).map((channel, index) => (
                    <TextField
                        key={channel}
                        size="small"
                        type="number"
                        label={channel}
                        value={value[index] ?? 0}
                        onChange={(event) => {
                            const next = [...value] as [number, number, number, number];
                            next[index] = Number(event.target.value);
                            onChange(next);
                        }}
                    />
                ))}
            </div>
        </div>
    );
}

export default function ShapesLayerPanel({
    config,
    updateLayer,
    association,
    availableFields,
    chartColorBy,
}: Props) {
    const fillColor = config.fillColor ?? [200, 200, 200, 120];
    const strokeColor = config.strokeColor ?? [255, 255, 255, 200];
    const tooltipFields = config.tooltipFields ?? [];
    const fillByColumn = config.fillColorByColumn?.columnName ?? "";

    const options = useMemo(() => {
        const set = new Set(availableFields);
        if (chartColorBy) set.add(chartColorBy);
        for (const field of tooltipFields) set.add(field);
        if (fillByColumn) set.add(fillByColumn);
        return [...set].sort();
    }, [availableFields, chartColorBy, fillByColumn, tooltipFields]);

    return (
        <div className="grid gap-3">
            {association.status === "resolved" && association.tableName && (
                <Typography variant="caption" color="text.secondary">
                    Associated table: {association.tableName}
                </Typography>
            )}
            {association.status === "ambiguous" && (
                <Typography variant="caption" color="warning.main">
                    Table association is ambiguous. Choose columns manually.
                </Typography>
            )}
            <div className="flex items-center gap-3">
                <span className="w-24 text-xs uppercase tracking-wide text-[hsl(var(--muted-foreground))]">
                    Stroke width
                </span>
                <Slider
                    size="small"
                    min={0}
                    max={8}
                    step={0.5}
                    value={config.strokeWidth ?? 1}
                    onChange={(_, value) => {
                        if (typeof value === "number") {
                            updateLayer({ strokeWidth: value });
                        }
                    }}
                />
            </div>
            <ColorFields
                label="Fill color"
                value={fillColor}
                onChange={(next) => updateLayer({ fillColor: next })}
            />
            <ColorFields
                label="Stroke color"
                value={strokeColor}
                onChange={(next) => updateLayer({ strokeColor: next })}
            />
            <Autocomplete
                multiple
                size="small"
                options={options}
                value={tooltipFields}
                onChange={(_, value) => updateLayer({ tooltipFields: value })}
                renderTags={(value, getTagProps) =>
                    value.map((option, index) => {
                        const { key, ...tagProps } = getTagProps({ index });
                        return <Chip key={key} {...tagProps} label={option} size="small" />;
                    })
                }
                renderInput={(params) => (
                    <TextField {...params} label="Tooltip fields" placeholder="Select columns" />
                )}
            />
            <Select
                size="small"
                displayEmpty
                value={fillByColumn}
                onChange={(event) => {
                    const columnName = event.target.value;
                    if (!columnName) {
                        updateLayer({ fillColorByColumn: undefined });
                        return;
                    }
                    updateLayer({
                        fillColorByColumn: {
                            columnName,
                            mode: "categorical",
                        },
                    });
                }}
            >
                <MenuItem value="">
                    <em>Static fill color</em>
                </MenuItem>
                {options.map((field) => (
                    <MenuItem key={field} value={field}>
                        {field}
                    </MenuItem>
                ))}
            </Select>
        </div>
    );
}
