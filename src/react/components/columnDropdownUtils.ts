import type { DataColumn, DataType } from "@/charts/charts";
import { columnMatchesType, type CTypes } from "@/lib/columnTypeHelpers";

export function getSelectableColumns(
    columns: DataColumn<DataType>[],
    type: CTypes | undefined,
    exclude?: string[],
) {
    return columns
        .filter((column) => !exclude?.includes(column.name))
        .filter((column) => !column.field.includes("|"))
        .filter((column) => columnMatchesType(column, type));
}
