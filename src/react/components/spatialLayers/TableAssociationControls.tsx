import { TextField, Typography } from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";
import type { ReactNode } from "react";

import type DataStore from "@/datastore/DataStore";
import type { FieldSpec, FieldSpecs } from "@/lib/columnTypeHelpers";
import { flattenFields } from "@/lib/columnTypeHelpers";
import type { TableAssociation } from "@/react/spatialdata/table_association";
import ColumnSelectionComponent from "../ColumnSelectionComponent";

type ShapesLayerConfig = Extract<LayerConfig, { type: "shapes" }>;

export type LayerFillColorByColumn = ShapesLayerConfig["fillColorByColumn"];

type Props = {
    association: TableAssociation;
    dataStore: DataStore;
    tooltipFields: string[];
    /** MDV's own record of the choice, when there is one — see {@link MdvFieldSpecs}. */
    tooltipFieldsSpec?: FieldSpecs;
    onTooltipFieldsChange: (next: string[], spec: FieldSpecs) => void;
};

type FillColorByColumnProps = {
    association: TableAssociation;
    dataStore: DataStore;
    fillColorByColumn?: LayerFillColorByColumn;
    fillColorByColumnSpec?: FieldSpec;
    onChange: (next: LayerFillColorByColumn | undefined, spec: FieldSpec | undefined) => void;
};

function ControlLabel({ children }: { children: string }) {
    return (
        <Typography
            fontSize="small"
            sx={{
                alignSelf: "center",
                justifySelf: "end",
                textAlign: "right",
                paddingRight: 2,
            }}
        >
            {children}
        </Typography>
    );
}

function LabeledControl({
    label,
    children,
}: {
    label: string;
    children: ReactNode;
}) {
    return (
        <div className="grid grid-cols-[7rem_minmax(0,1fr)] items-center gap-2">
            <ControlLabel>{label}</ControlLabel>
            <div className="min-w-0">{children}</div>
        </div>
    );
}

export default function TableAssociationControls({
    association,
    dataStore,
    tooltipFields,
    tooltipFieldsSpec,
    onTooltipFieldsChange,
}: Props) {
    return (
        <div className="grid gap-2">
            {association.status === "resolved" && association.tableName && (
                <Typography variant="caption" color="text.secondary">
                    Associated datasource: {association.dataSourceName} / {association.tableName}
                    {association.featureCount !== undefined &&
                        association.matchedFeatureCount !== undefined &&
                        ` (${association.matchedFeatureCount}/${association.featureCount} features)`}
                </Typography>
            )}
            {association.status === "loading" && (
                <Typography variant="caption" color="text.secondary">
                    Resolving table association...
                </Typography>
            )}
            {association.status === "ambiguous" && (
                <Typography variant="caption" color="warning.main">
                    Table association is ambiguous. Choose columns manually.
                </Typography>
            )}
            {association.status === "none" && (
                <Typography variant="caption" color="text.secondary">
                    No associated table inferred.
                </Typography>
            )}
            <LabeledControl label="Tooltip fields">
                <ColumnSelectionComponent
                    multiple
                    dataStore={dataStore}
                    // The spec when there is one: it is what the picker gave us, so it
                    // is what the picker can read back — an active link shown as an
                    // active link rather than as the columns it happens to resolve to.
                    current_value={tooltipFieldsSpec ?? tooltipFields}
                    placeholder="Select columns"
                    setSelectedColumn={(next) => onTooltipFieldsChange(flattenFields(next), next)}
                />
            </LabeledControl>
        </div>
    );
}

/**
 * What we claim about a colour column when writing it into the render stack:
 * nothing. `"auto"` lets the reader decide categorical-vs-continuous from the
 * column's declared kind, which is the only party that knows it reliably.
 *
 * This was `"categorical"` for every column, which is how a `double` ended up
 * with one palette entry per distinct float instead of a ramp. An explicit mode
 * is honoured verbatim by the viewer — it exists so a user can force categorical
 * on, say, integer cluster codes — so asserting one here overrode the detection
 * rather than informing it. If a per-layer override is ever wanted, it belongs in
 * the UI as a choice, not as a constant.
 *
 * For a column in the table's obs this is overwritten before the config reaches
 * the viewer: `fillColorSchemeFromDataStore` sets the mode alongside the palette
 * it derived, because a `byValue` palette is meaningless under a continuous mode.
 * What is written here is what a column OUTSIDE obs carries, and what a reader of
 * the saved stack falls back to.
 */
const FILL_COLOR_MODE = "auto" as const;

export function FillColorByColumnControl({
    association,
    dataStore,
    fillColorByColumn,
    fillColorByColumnSpec,
    onChange,
}: FillColorByColumnProps) {
    const fillByColumn = fillColorByColumn?.columnName;

    return (
        <div className="grid gap-2">
            <LabeledControl label="Fill color">
                {association.status === "resolved" ? (
                    <ColumnSelectionComponent
                        multiple={false}
                        dataStore={dataStore}
                        current_value={fillColorByColumnSpec ?? fillByColumn}
                        placeholder="Select column"
                        optional
                        clearSelectedColumn={() => onChange(undefined, undefined)}
                        setSelectedColumn={(selection) => {
                            // Not necessarily a column name: the picker's "active link"
                            // tab returns a query, whose column is whatever the linked
                            // datasource has selected right now. `flattenFields` reads
                            // both, and an empty result means the link has not finished
                            // initialising — the spec is still stored, and the projection
                            // fills the colour column in when it resolves.
                            const columnName = flattenFields(selection)[0];
                            onChange(
                                columnName ? { columnName, mode: FILL_COLOR_MODE } : undefined,
                                selection,
                            );
                        }}
                    />
                ) : (
                    <TextField
                        size="small"
                        disabled
                        value="Static fill color"
                        slotProps={{ input: { readOnly: true } }}
                    />
                )}
            </LabeledControl>
        </div>
    );
}
