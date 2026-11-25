import { Typography } from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";
import { observer } from "mobx-react-lite";
import type { ReactNode } from "react";
import type { GridBySpec } from "@/charts/charts";

import type DataStore from "@/datastore/DataStore";
import type { FieldSpec, FieldSpecs } from "@/lib/columnTypeHelpers";
import { flattenFields } from "@/lib/columnTypeHelpers";
import type { TableAssociation } from "@/react/spatialdata/table_association";
import ColumnSelectionComponent from "../ColumnSelectionComponent";
import GridByControl from "../GridByControl";

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
    fillColorByColumnSpec?: FieldSpec | FieldSpecs;
    onChange: (next: LayerFillColorByColumn | undefined, spec: FieldSpec | FieldSpecs | undefined) => void;
    /** The chart's grid by this layer's fill colour. */
    gridBy?: GridBySpec;
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

function LabeledControl({ label, children }: { label: string; children: ReactNode }) {
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
                    More than one table annotates this element, so there are no columns to show or colour by.
                </Typography>
            )}
            {association.status === "none" && (
                <Typography variant="caption" color="text.secondary">
                    No table annotates this element: tooltips show the feature id only.
                </Typography>
            )}
            {/* Fields come from the annotating table's datasource; without one there is nothing to read. */}
            {association.status === "resolved" && (
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
            )}
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

/**
 * The layer's colour column. With `gridBy`, a control to make the chart a grid of it: of a cell
 * per category of a categorical column, or of a cell per numeric column, when the picker takes
 * several (the layer itself takes the first). Turned off, it shows the first again.
 */
export const FillColorByColumnControl = observer(function FillColorByColumnControl({
    association,
    dataStore,
    fillColorByColumn,
    fillColorByColumnSpec,
    onChange,
    gridBy,
}: FillColorByColumnProps) {
    const fillByColumn = fillColorByColumn?.columnName;
    const split = gridBy?.get() ?? null;
    const available = gridBy?.available() ?? null;

    // Without a table there is no column to colour by: the layer's own colour controls apply.
    if (association.status !== "resolved") return null;
    // Not necessarily column names: the picker's "active link" tab returns a query, whose
    // columns are whatever the linked datasource has selected right now. `flattenFields` reads
    // both, and an empty result means the link has not finished initialising — the spec is
    // still stored, and the projection fills the colour column in when it resolves.
    const select = (selection: FieldSpec | FieldSpecs) => {
        const columnName = flattenFields(selection)[0];
        onChange(columnName ? { columnName, mode: FILL_COLOR_MODE } : undefined, selection);
        // the split follows the new column's type
        if (gridBy?.get()) gridBy.set(true);
    };
    const specs = Array.isArray(fillColorByColumnSpec)
        ? fillColorByColumnSpec
        : fillColorByColumnSpec
          ? [fillColorByColumnSpec]
          : fillByColumn
            ? [fillByColumn]
            : [];
    return (
        <div className="grid gap-2">
            <LabeledControl label="Colour by">
                <div className="grid gap-1">
                    {split === "columns" ? (
                        <ColumnSelectionComponent
                            multiple
                            type="number"
                            dataStore={dataStore}
                            current_value={specs}
                            placeholder="Select columns"
                            setSelectedColumn={select}
                        />
                    ) : (
                        <ColumnSelectionComponent
                            multiple={false}
                            dataStore={dataStore}
                            current_value={specs[0]}
                            placeholder="Select column"
                            optional
                            clearSelectedColumn={() => onChange(undefined, undefined)}
                            setSelectedColumn={select}
                        />
                    )}
                    {gridBy && <GridByControl split={split} available={available} onChange={gridBy.set} />}
                </div>
            </LabeledControl>
        </div>
    );
});
