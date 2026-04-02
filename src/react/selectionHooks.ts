import { useEffect, useState, useCallback, useId } from "react";
import { useChart, useDataStore } from "./context";

export type HighlightModifiers = {
    add?: boolean;
    remove?: boolean;
    toggle?: boolean;
};

export type ModifierEventLike = {
    shiftKey?: boolean;
    ctrlKey?: boolean;
    metaKey?: boolean;
    srcEvent?: unknown;
    sourceEvent?: unknown;
} | null | undefined;

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
    // Create a copy to avoid a reference to mutable state - not that we anticipate an actual issue in this case.
    const [highlightedIndices, setHighlightedIndices] = useState<number[]>(
        dataStore.highightedData?.slice() || []
    );
    const listenerId = useId();

    useEffect(() => {
        const listener = (type: string, data: { indexes: number[] | Record<number, number> }) => {
            if (type === "data_highlighted") {
                // Handle both array and object formats (is this necessary?)
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

export function getHighlightModifiers(event?: ModifierEventLike): HighlightModifiers {
    const raw = event ?? {};
    const nested =
        typeof raw === "object" && raw !== null
            ? ((raw.srcEvent ?? raw.sourceEvent) as ModifierEventLike)
            : undefined;
    const resolved = nested ?? raw;
    return {
        add: !!resolved?.shiftKey,
        toggle: !!(resolved?.ctrlKey || resolved?.metaKey),
    };
}

export function useHighlightRowUpdater() {
    const highlightedIndices = useHighlightedIndices();
    const highlightRows = useHighlightRows();

    return useCallback(
        (rowIndex: number, modifiers?: HighlightModifiers) => {
            if (rowIndex < 0) return;
            const { add = false, remove = false, toggle = false } = modifiers ?? {};
            if (toggle) {
                const nextHighlights = new Set(highlightedIndices);
                if (nextHighlights.has(rowIndex)) nextHighlights.delete(rowIndex);
                else nextHighlights.add(rowIndex);
                highlightRows(Array.from(nextHighlights));
                return;
            }
            if (remove) {
                if (!highlightedIndices.includes(rowIndex)) return;
                highlightRows(highlightedIndices.filter((index) => index !== rowIndex));
                return;
            }
            if (add) {
                if (highlightedIndices.includes(rowIndex)) return;
                highlightRows([...highlightedIndices, rowIndex]);
                return;
            }
            highlightRows([rowIndex]);
        },
        [highlightRows, highlightedIndices],
    );
}
