import { useEffect, useState, useCallback, useId } from "react";
import { useChart, useDataStore } from "./context";

/**
 * Hook to get the index of the highlighted data point.
 * Delegates to useHighlightedIndices() and returns the first index.
 */
export function useHighlightedIndex() {
    const highlightedIndices = useHighlightedIndices();
    return highlightedIndices[0] ?? -1;
}

/**
 * Hook to get all highlighted indices as an array.
 * Supports multiple highlighted rows.
 */
export function useHighlightedIndices() {
    const dataStore = useDataStore();
    const [highlightedIndices, setHighlightedIndices] = useState<number[]>(dataStore.highightedData || []);
    const listenerId = useId();

    useEffect(() => {
        const listener = (type: string, data: { indexes: number[] | { [k: number]: number } }) => {
            if (type === "data_highlighted") {
                // Handle both array and object formats
                const indices = Array.isArray(data.indexes) 
                    ? data.indexes 
                    : Object.values(data.indexes);
                setHighlightedIndices(indices);
            }
        };

        dataStore.addListener(listenerId, listener);

        return () => {
            dataStore.removeListener(listenerId);
        };
    }, [dataStore, listenerId]);

    return highlightedIndices;
}

/**
 * Hook that returns a function to highlight rows.
 * The returned function calls dataStore.dataHighlighted with the provided row indices.
 */
export function useHighlightRows() {
    const chart = useChart();
    const dataStore = useDataStore();
    
    return useCallback((rowIndexes: number[]) => {
        // dataHighlighted expects an array of numbers
        dataStore.dataHighlighted(rowIndexes, chart);
    }, [chart, dataStore]);
}
