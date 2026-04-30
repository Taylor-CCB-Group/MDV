import { useEffect, useState, useCallback, useRef } from "react";
import { useChart, useDataStore } from "./context";

type CategoryRows = {
    rows: ArrayLike<number>;
};

let nextHighlightedIndicesListenerId = 0;

export type HighlightModifierInput = {
    shiftKey?: boolean;
    ctrlKey?: boolean;
    metaKey?: boolean;
};

export type HighlightModifierState = {
    add: boolean;
    toggle: boolean;
};

function getUniqueRows(rowIndexes: ArrayLike<number>) {
    return Array.from(new Set(Array.from(rowIndexes)));
}

/**
 * Normalizes keyboard modifiers into the shared highlight semantics used by row-selection UIs.
 * Ctrl/Cmd toggle wins over Shift-add when both are pressed.
 */
export function getHighlightModifierState(
    modifiers?: HighlightModifierInput,
): HighlightModifierState {
    const toggle = !!(modifiers?.ctrlKey || modifiers?.metaKey);
    return {
        add: !toggle && !!modifiers?.shiftKey,
        toggle,
    };
}

export function getNextHighlightedRows(
    rowIndexes: ArrayLike<number>,
    currentHighlightedRows: ArrayLike<number>,
    modifiers?: HighlightModifierInput,
) {
    const targetRows = getUniqueRows(rowIndexes);
    const { add, toggle } = getHighlightModifierState(modifiers);

    if (toggle) {
        const currentRows = getUniqueRows(currentHighlightedRows);
        const currentRowSet = new Set(currentRows);
        const shouldRemoveGroup = targetRows.every((rowIndex) =>
            currentRowSet.has(rowIndex),
        );

        if (shouldRemoveGroup) {
            const targetRowSet = new Set(targetRows);
            return currentRows.filter((rowIndex) => !targetRowSet.has(rowIndex));
        }

        return getUniqueRows([...currentRows, ...targetRows]);
    }

    if (add) {
        return getUniqueRows([...Array.from(currentHighlightedRows), ...targetRows]);
    }

    return targetRows;
}

export function getHighlightedRowsFromCategoryIndices(
    categoryIndices: number[],
    categories: CategoryRows[],
) {
    const highlightedRows = new Set<number>();
    for (const categoryIndex of categoryIndices) {
        const categoryRows = categories[categoryIndex]?.rows;
        if (!categoryRows) continue;
        for (const row of Array.from(categoryRows)) {
            highlightedRows.add(row);
        }
    }
    return Array.from(highlightedRows);
}

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
    const listenerKeyRef = useRef<string | null>(null);
    let listenerKey = listenerKeyRef.current;
    if (listenerKey === null) {
        listenerKey = `highlighted-indices-${nextHighlightedIndicesListenerId++}`;
        listenerKeyRef.current = listenerKey;
    }

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

        dataStore.addListener(listenerKey, listener);

        return () => {
            dataStore.removeListener(listenerKey);
        };
    }, [dataStore, listenerKey]);

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

/**
 * Hook that highlights rows while respecting the shared replace/add/toggle keyboard modifiers.
 */
export function useHighlightRowsWithModifiers() {
    const highlightRows = useHighlightRows();
    const highlightedIndices = useHighlightedIndices();

    return useCallback(
        (rowIndexes: ArrayLike<number>, modifiers?: HighlightModifierInput) => {
            highlightRows(
                getNextHighlightedRows(rowIndexes, highlightedIndices, modifiers),
            );
        },
        [highlightRows, highlightedIndices],
    );
}

/**
 * Hook that highlights rows by selected category indices.
 * This is useful when the interaction model is category-based but the shared
 * application highlight state is still expressed in row indices.
 */
export function useHighlightCategoryRows(categories: CategoryRows[]) {
    const highlightRows = useHighlightRows();

    return useCallback((categoryIndices: number[]) => {
        const rows = getHighlightedRowsFromCategoryIndices(categoryIndices, categories);
        highlightRows(rows);
    }, [categories, highlightRows]);
}
