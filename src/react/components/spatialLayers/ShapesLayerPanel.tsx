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

const toHex = (value: [number, number, number, number]) =>
    `#${value
        .slice(0, 3)
        .map((channel) => Math.max(0, Math.min(255, channel)).toString(16).padStart(2, "0"))
        .join("")}`;

const fromHex = (hex: string): [number, number, number] => [
    Number.parseInt(hex.slice(1, 3), 16),
    Number.parseInt(hex.slice(3, 5), 16),
    Number.parseInt(hex.slice(5, 7), 16),
];

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
            <input
                type="color"
                value={toHex(value)}
                onChange={(event) => {
                    const [r, g, b] = fromHex(event.target.value);
                    onChange([r, g, b, value[3] ?? 255]);
                }}
            />
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
                    //nb seems like this is being capped by a default config.strokeWidthMaxPixels = 1
                    //we should change that (could overwrite it ourselves) but for now limit this range
                    max={1}
                    step={0.1}
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
