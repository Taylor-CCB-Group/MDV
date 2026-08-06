import { TextField, Typography } from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";
import type { ReactNode } from "react";

import type DataStore from "@/datastore/DataStore";
import type { FieldSpecs } from "@/lib/columnTypeHelpers";
import type { TableAssociation } from "@/react/spatialdata/table_association";
import ColumnSelectionComponent from "../ColumnSelectionComponent";

type ShapesLayerConfig = Extract<LayerConfig, { type: "shapes" }>;

export type LayerFillColorByColumn = ShapesLayerConfig["fillColorByColumn"];

type Props = {
    association: TableAssociation;
    dataStore: DataStore;
    tooltipFields: string[];
    onTooltipFieldsChange: (next: string[]) => void;
};

type FillColorByColumnProps = {
    association: TableAssociation;
    dataStore: DataStore;
    fillColorByColumn?: LayerFillColorByColumn;
    onChange: (next: LayerFillColorByColumn | undefined) => void;
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

function isConcreteFieldName(column: FieldSpecs[number]): column is string {
    return typeof column === "string";
}

function concreteFieldNames(columns: FieldSpecs): string[] {
    return columns.filter(isConcreteFieldName);
}

export default function TableAssociationControls({
    association,
    dataStore,
    tooltipFields,
    onTooltipFieldsChange,
}: Props) {
    return (
        <div className="grid gap-2">
            {association.status === "resolved" && association.tableName && (
                <Typography variant="caption" color="text.secondary">
                    Associated datasource: {association.dataSourceName} /{" "}
                    {association.tableName}
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
                    current_value={tooltipFields}
                    placeholder="Select columns"
                    setSelectedColumn={(next) => onTooltipFieldsChange(concreteFieldNames(next))}
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
 * MDV renders these colours itself (see `withoutViewerFillColorColumn`), keyed on
 * the DataStore column's own datatype, so this value is about what a *reader* of
 * the saved stack sees — not about what MDV draws today.
 */
const FILL_COLOR_MODE = "auto" as const;

export function FillColorByColumnControl({
    association,
    dataStore,
    fillColorByColumn,
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
                        current_value={fillByColumn}
                        placeholder="Select column"
                        optional
                        clearSelectedColumn={() => onChange(undefined)}
                        setSelectedColumn={(columnName) => {
                            if (typeof columnName !== "string") return;
                            onChange({
                                columnName,
                                mode: FILL_COLOR_MODE,
                            });
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
