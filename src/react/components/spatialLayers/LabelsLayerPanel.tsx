import {
    Autocomplete,
    Chip,
    TextField,
    Typography,
} from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";

type LabelsLayerConfig = Extract<LayerConfig, { type: "labels" }>;
import { useMemo } from "react";
import type { TableAssociation } from "@/react/spatialdata/table_association";

type Props = {
    config: LabelsLayerConfig;
    updateLayer: (updates: Partial<LabelsLayerConfig>) => void;
    association: TableAssociation;
    availableFields: string[];
};

export default function LabelsLayerPanel({
    config,
    updateLayer,
    association,
    availableFields,
}: Props) {
    const selected = config.tooltipFields ?? [];
    const options = useMemo(() => {
        const set = new Set(availableFields);
        for (const field of selected) set.add(field);
        return [...set].sort();
    }, [availableFields, selected]);

    return (
        <div className="grid gap-2">
            {association.status === "ambiguous" && (
                <Typography variant="caption" color="warning.main">
                    Multiple tables match this element. Tooltip fields may be incomplete.
                </Typography>
            )}
            {association.status === "none" && (
                <Typography variant="caption" color="text.secondary">
                    No associated table inferred. Add tooltip fields manually.
                </Typography>
            )}
            <Autocomplete
                multiple
                size="small"
                options={options}
                value={selected}
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
        </div>
    );
}
