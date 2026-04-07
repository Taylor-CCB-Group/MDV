import { useDataStore } from "@/react/context";
import { useEffect, useId, useMemo, useState } from "react";
import {
    getAvailableCategorySelectionValues,
    shouldRefreshCategorySelectionValues,
} from "./categorySelectionUtils";

/**
 * Returns the current selectable category values for a source column and refreshes
 * them when datastore events indicate that the column's available categories may
 * have changed in place.
 */
export function useLiveCategorySelectionValues(
    sourceColumn: string | undefined,
) {
    const dataStore = useDataStore();
    const listenerId = useId();
    const [refreshVersion, setRefreshVersion] = useState(0);

    useEffect(() => {
        const key = `CategorySelectionOptions-${listenerId}`;
        const listener = (
            eventType: string,
            eventData: { columns?: string[] } | string | number | undefined,
        ) => {
            if (
                !shouldRefreshCategorySelectionValues(
                    eventType,
                    eventData,
                    sourceColumn,
                )
            ) {
                return;
            }
            setRefreshVersion((version) => version + 1);
        };

        dataStore.addListener(key, listener);
        return () => dataStore.removeListener(key);
    }, [dataStore, listenerId, sourceColumn]);

    return useMemo(
        () => getAvailableCategorySelectionValues(dataStore, sourceColumn),
        [dataStore, refreshVersion, sourceColumn],
    );
}
