import TableChartOutlinedIcon from "@mui/icons-material/TableChartOutlined";
import { Typography } from "@mui/material";
import type { RenderStackSpatialElementType } from "@spatialdata/layers";
import { useDataSources } from "@/react/hooks";
import { useChartSpatialData } from "@/react/spatialdata/shared_spatial_data";
import {
    type DataSourceAssociationCandidate,
    resolveAssociatedElementTable,
} from "@/react/spatialdata/table_association";

/** Which table annotates an element, as the layer list and add menu show it. */
export type ElementTableSummary =
    | { status: "none" }
    | { status: "ambiguous"; detail: string }
    | { status: "resolved"; tableName: string; dataSourceName: string };

type AssociationSource = Parameters<typeof resolveAssociatedElementTable>[0]["spatialData"];

/**
 * Whether a table annotates the element, by the same association the layer panels use
 * (`resolveAssociatedElementTable`), without loading the element: it reads the store's
 * metadata and the datasources' provenance only, so it is cheap enough for every row.
 *
 * Only shapes and labels: they are the layers MDV colours by a column and gives tooltip
 * fields. A table on a points or image element is not used by their panels, so it is not
 * advertised: those get `null`, as does any element while the store is being read.
 */
export function summarizeElementTable<TDataStore>(
    spatialData: AssociationSource,
    dataSources: DataSourceAssociationCandidate<TDataStore>[],
    elementType: RenderStackSpatialElementType,
    elementKey: string,
): ElementTableSummary | null {
    if (elementType !== "shapes" && elementType !== "labels") return null;
    if (!spatialData) return null;
    const table = resolveAssociatedElementTable({ spatialData, elementType, elementKey, dataSources });
    switch (table.status) {
        case "none":
            return { status: "none" };
        case "ambiguous":
            return {
                status: "ambiguous",
                detail: table.dataSourceNames
                    ? `table '${table.tableNames[0]}' matches several datasources (${table.dataSourceNames.join(", ")})`
                    : `several tables annotate it (${table.tableNames.join(", ")})`,
            };
        case "resolved":
            return { status: "resolved", tableName: table.tableName, dataSourceName: table.dataSourceName };
    }
}

/** {@link summarizeElementTable} for the chart's store and the project's datasources. */
export function useElementTableSummary() {
    const { spatialData } = useChartSpatialData();
    const dataSources = useDataSources();
    return (elementType: RenderStackSpatialElementType, elementKey: string) =>
        summarizeElementTable(spatialData ?? undefined, dataSources, elementType, elementKey);
}

/** What a summary means, for a tooltip. */
export function elementTableDescription(summary: ElementTableSummary): string {
    switch (summary.status) {
        case "none":
            return "No table annotates this element: it cannot be coloured by a column, and tooltips show the feature id only.";
        case "ambiguous":
            return `Ambiguous table: ${summary.detail}, so there are no columns to colour by or show.`;
        case "resolved":
            return `Annotated by datasource '${summary.dataSourceName}' (table '${summary.tableName}'): colour by its columns and show them in tooltips.`;
    }
}

/** A table icon for an annotated element (amber when the association is ambiguous); nothing otherwise. */
export function ElementTableIcon({ summary }: { summary: ElementTableSummary }) {
    if (summary.status === "none") return null;
    const description = elementTableDescription(summary);
    return (
        <span className="flex flex-shrink-0" title={description} role="img" aria-label={description}>
            <TableChartOutlinedIcon
                fontSize="small"
                color={summary.status === "ambiguous" ? "warning" : "action"}
                sx={{ fontSize: 16 }}
            />
        </span>
    );
}

/** The add menu's note on an element: the annotating datasource, "ambiguous table", or "no table". */
export function ElementTableMenuHint({ summary }: { summary: ElementTableSummary }) {
    return (
        <span className="ml-auto pl-4 flex items-center gap-1" title={elementTableDescription(summary)}>
            <ElementTableIcon summary={summary} />
            <Typography
                variant="caption"
                color={summary.status === "ambiguous" ? "warning.main" : "text.secondary"}
                sx={{ fontStyle: summary.status === "none" ? "italic" : undefined }}
            >
                {summary.status === "resolved"
                    ? summary.dataSourceName
                    : summary.status === "ambiguous"
                      ? "ambiguous table"
                      : "no table"}
            </Typography>
        </span>
    );
}
