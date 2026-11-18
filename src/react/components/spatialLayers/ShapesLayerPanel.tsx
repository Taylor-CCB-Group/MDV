import { Slider } from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";

import type DataStore from "@/datastore/DataStore";
import { mdvFieldSpecsOf, type WithMdvFieldSpecs } from "@/react/spatialdata/field_spec_projection";
import type { TableAssociation } from "@/react/spatialdata/table_association";
import { ColorFields } from "./ColorFields";
import TableAssociationControls, { FillColorByColumnControl } from "./TableAssociationControls";

type ShapesLayerConfig = WithMdvFieldSpecs<Extract<LayerConfig, { type: "shapes" }>>;

type Props = {
    config: ShapesLayerConfig;
    updateLayer: (updates: Partial<ShapesLayerConfig>) => void;
    association: TableAssociation;
    dataStore: DataStore;
};

export default function ShapesLayerPanel({ config, updateLayer, association, dataStore }: Props) {
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
                label="Fill colour"
                rgb={[fillColor[0], fillColor[1], fillColor[2]]}
                alpha={fillColor[3] ?? 255}
                caption={
                    config.fillColorByColumn
                        ? "Coloured by column; the opacity here still applies."
                        : association.status === "resolved"
                          ? "Choose a column above to colour by its values instead."
                          : undefined
                }
                onChange={([r, g, b], a) => updateLayer({ fillColor: [r, g, b, a ?? 255] })}
            />
            <ColorFields
                label="Stroke colour"
                rgb={[strokeColor[0], strokeColor[1], strokeColor[2]]}
                alpha={strokeColor[3] ?? 255}
                onChange={([r, g, b], a) => updateLayer({ strokeColor: [r, g, b, a ?? 255] })}
            />
        </div>
    );
}
