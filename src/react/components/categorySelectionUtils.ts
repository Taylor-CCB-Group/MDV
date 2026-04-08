import type DataStore from "@/datastore/DataStore";
import { getMultitextItemValues, inferMultitextEmptyItem } from "@/lib/multitext";

export function getAvailableCategorySelectionValues(
    dataStore: Pick<DataStore, "columnIndex" | "getColumnValues">,
    sourceColumn: string | undefined,
) {
    if (!sourceColumn) {
        return [];
    }

    try {
        const source = dataStore.columnIndex[sourceColumn];
        if (source?.datatype === "multitext") {
            return getMultitextItemValues(source, {
                emptyItem: inferMultitextEmptyItem(source),
            }).slice();
        }
        return dataStore.getColumnValues(sourceColumn).slice();
    } catch (error) {
        console.error(
            `error updating category selection values with '${sourceColumn}' (${error})`,
        );
        return [];
    }
}

export function getCategorySelectionDropdownKey(
    sourceColumn: string | undefined,
    availableValues: string[],
) {
    return `${sourceColumn || ""}\u0000${availableValues.join("\u0000")}`;
}
