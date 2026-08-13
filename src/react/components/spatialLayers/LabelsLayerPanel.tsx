import type { LayerConfig } from "@spatialdata/vis";

import type DataStore from "@/datastore/DataStore";
import { mdvFieldSpecsOf, type WithMdvFieldSpecs } from "@/react/spatialdata/field_spec_projection";
import type { TableAssociation } from "@/react/spatialdata/table_association";
import TableAssociationControls, {
    FillColorByColumnControl,
    type LayerFillColorByColumn,
} from "./TableAssociationControls";

type LabelsLayerConfig = WithMdvFieldSpecs<
    Extract<LayerConfig, { type: "labels" }> & {
        fillColorByColumn?: LayerFillColorByColumn;
    }
>;

type Props = {
    config: LabelsLayerConfig;
    updateLayer: (updates: Partial<LabelsLayerConfig>) => void;
    association: TableAssociation;
    dataStore: DataStore;
};

export default function LabelsLayerPanel({
    config,
    updateLayer,
    association,
    dataStore,
}: Props) {
    const specs = mdvFieldSpecsOf(config);

    return (
        <div className="grid gap-2">
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
        </div>
    );
}
