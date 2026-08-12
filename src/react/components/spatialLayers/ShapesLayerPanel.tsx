import { Slider, Typography } from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";

import type DataStore from "@/datastore/DataStore";
import { mdvFieldSpecsOf, type WithMdvFieldSpecs } from "@/react/spatialdata/field_spec_projection";
import type { TableAssociation } from "@/react/spatialdata/table_association";
import TableAssociationControls, { FillColorByColumnControl } from "./TableAssociationControls";

type ShapesLayerConfig = WithMdvFieldSpecs<Extract<LayerConfig, { type: "shapes" }>>;

type Props = {
    config: ShapesLayerConfig;
    updateLayer: (updates: Partial<ShapesLayerConfig>) => void;
    association: TableAssociation;
    dataStore: DataStore;
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
            <div className="flex items-center gap-3">
                <input
                    type="color"
                    value={toHex(value)}
                    onChange={(event) => {
                        const [r, g, b] = fromHex(event.target.value);
                        onChange([r, g, b, value[3] ?? 255]);
                    }}
                />
                <Slider
                    size="small"
                    min={0}
                    max={255}
                    step={1}
                    value={value[3] ?? 255}
                    onChange={(_, next) => {
                        if (typeof next === "number") {
                            onChange([value[0], value[1], value[2], next]);
                        }
                    }}
                />
            </div>
        </div>
    );
}

export default function ShapesLayerPanel({
    config,
    updateLayer,
    association,
    dataStore,
}: Props) {
    const fillColor = config.fillColor ?? [200, 200, 200, 120];
    const strokeColor = config.strokeColor ?? [255, 255, 255, 200];
    const tooltipFields = config.tooltipFields ?? [];
    const specs = mdvFieldSpecsOf(config);

    return (
        <div className="grid gap-3">
            <TableAssociationControls
                association={association}
                dataStore={dataStore}
                tooltipFields={tooltipFields}
                tooltipFieldsSpec={specs?.tooltipFields}
                onTooltipFieldsChange={(next, spec) =>
                    updateLayer({
                        tooltipFields: next,
                        mdvFieldSpecs:{ ...specs, tooltipFields: spec },
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
                        mdvFieldSpecs:{ ...specs, fillColorByColumn: spec },
                    })
                }
            />
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
        </div>
    );
}
