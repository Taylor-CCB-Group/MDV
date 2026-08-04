import type { LayerConfig } from "@spatialdata/vis";

import type DataStore from "@/datastore/DataStore";
import type { TableAssociation } from "@/react/spatialdata/table_association";
import TableAssociationControls, {
    FillColorByColumnControl,
    type LayerFillColorByColumn,
} from "./TableAssociationControls";

type LabelsLayerConfig = Extract<LayerConfig, { type: "labels" }> & {
    fillColorByColumn?: LayerFillColorByColumn;
};

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
    return (
        <div className="grid gap-2">
            <TableAssociationControls
                association={association}
                dataStore={dataStore}
                tooltipFields={config.tooltipFields ?? []}
                onTooltipFieldsChange={(next) => updateLayer({ tooltipFields: next })}
            />
            <FillColorByColumnControl
                association={association}
                dataStore={dataStore}
                fillColorByColumn={config.fillColorByColumn}
                onChange={(next) => updateLayer({ fillColorByColumn: next })}
            />
        </div>
    );
}
