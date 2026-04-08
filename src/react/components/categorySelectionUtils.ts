import type DataStore from "@/datastore/DataStore";
import { getMultitextItemValues, inferMultitextEmptyItem } from "@/lib/multitext";

type CategorySelectionDataChange = {
    columns?: string[];
};

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

export function shouldRefreshCategorySelectionValues(
    eventType: string,
    eventData: CategorySelectionDataChange | string | number | undefined,
    sourceColumn: string | undefined,
) {
    if (!sourceColumn) {
        return false;
    }

    if (eventType === "data_added") {
        return true;
    }

    if (eventType === "column_removed") {
        return eventData === sourceColumn;
    }

    if (eventType !== "data_changed") {
        return false;
    }

    return typeof eventData === "object"
        && eventData !== null
        && Array.isArray(eventData.columns)
        && eventData.columns.includes(sourceColumn);
}
