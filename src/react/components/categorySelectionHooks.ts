import { useDataStore } from "@/react/context";
import { useEffect, useRef, useMemo, useState } from "react";
import {
    getAvailableCategorySelectionValues,
    shouldRefreshCategorySelectionValues,
} from "./categorySelectionUtils";

let nextCategorySelectionListenerId = 0;

/**
 * Returns the current selectable category values for a source column and refreshes
 * them when datastore events indicate that the column's available categories may
 * have changed in place.
 */
export function useLiveCategorySelectionValues(
    sourceColumn: string | undefined,
) {
    const dataStore = useDataStore();
    // const listenerId = useId(); // theoretical id collision with different react roots
    //https://github.com/Taylor-CCB-Group/MDV/pull/401#discussion_r3046010093
    const listenerKeyRef = useRef(
        `CategorySelectionOptions-${nextCategorySelectionListenerId++}`
    );
    const [refreshVersion, setRefreshVersion] = useState(0);

    useEffect(() => {
        const key = listenerKeyRef.current;
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
    }, [dataStore, listenerKeyRef, sourceColumn]);

    return useMemo(
        () => getAvailableCategorySelectionValues(dataStore, sourceColumn),
        [dataStore, refreshVersion, sourceColumn],
    );
}
