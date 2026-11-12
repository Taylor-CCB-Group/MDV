import type { DataColumn, DataType } from "@/charts/charts";
import { type CTypes, columnMatchesType } from "@/lib/columnTypeHelpers";

export function getSelectableColumns(columns: DataColumn<DataType>[], type: CTypes | undefined, exclude?: string[]) {
    return columns
        .filter((column) => !exclude?.includes(column.name))
        .filter((column) => !column.field.includes("|"))
        .filter((column) => columnMatchesType(column, type));
}
