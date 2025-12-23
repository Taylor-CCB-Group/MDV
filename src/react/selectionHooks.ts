import { useEffect, useState, useCallback } from "react";
import { useChart, useDataStore } from "./context";

/**
 * Hook to get the index of the highlighted data point.
 * May want to change how this works in terms of context etc, current implementation is simple
 * will not behave well if used multiple times by the same chart, etc.
 */
export function useHighlightedIndex() {
    const chart = useChart();
    // const { highlightedData } = chart.dataStore;
    const [highlightedIndex, setHighlightedIndex] = useState<number>(-1);

    useEffect(() => {
        console.log("useHighlightedIndex effect");
        // todo check the type - but it's definitely indexable by number, returning number
        chart.onDataHighlighted = ({ indexes }: { indexes: { [k: number]: number } }) => {
            setHighlightedIndex(indexes[0]);
        };
    }, [chart]);

    return highlightedIndex;
}

/**
 * Hook to get all highlighted indices as an array.
 * Supports multiple highlighted rows.
 */
export function useHighlightedIndices(): number[] {
    const chart = useChart();
    const [highlightedIndices, setHighlightedIndices] = useState<number[]>([]);

    useEffect(() => {
        chart.onDataHighlighted = ({ indexes }: { indexes: number[] | { [k: number]: number } }) => {
            // Handle both array and object formats
            const indices = Array.isArray(indexes) ? indexes : Object.values(indexes);
            setHighlightedIndices(indices);
        };
    }, [chart]);

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
