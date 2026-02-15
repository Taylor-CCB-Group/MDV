import type { DataType, LoadedDataColumn } from "@/charts/charts";
import { useCallback, useRef } from "react";
import type { OnBeforeEditCellEventArgs, OnCellChangeEventArgs, SlickgridReactInstance } from "slickgrid-react";
import { setCellValueFromString } from "../utils/valueReplacementUtil";
import type DataStore from "@/datastore/DataStore";
import type { FeedbackAlert } from "../components/FeedbackAlertComponent";

/**
 *
 * Hook for handling cell editing in the grid
 *
 * Uses refs, not values, for orderedParamColumns and sortedFilteredIndices
 * because SlickGrid stores callback references internally. Using refs
 * ensures the callbacks always access the current data avoiding stale closures.
 * 
 * - Stores the old value for comparison when applying edit
 * - Uses the util function setCellValueFromString for editing
 * - Calls dataStore.dataChanged for dimension filter to refilter
 * - Rerender the grid after updation
 * 
*/
const useEditCell = (
    orderedParamColumnsRef: React.MutableRefObject<LoadedDataColumn<DataType>[]>,
    sortedFilteredIndicesRef: React.MutableRefObject<Uint32Array>,
    dataStore: DataStore,
    gridRef: React.MutableRefObject<SlickgridReactInstance | null>,
    setFeedbackAlert: (alert: FeedbackAlert) => void,
) => {
    const oldCellValueRef = useRef<string | null>(null);

    // Store the old value of the cell in the oldCellValueRef
    const handleBeforeEditCell = useCallback(
        (
            e: CustomEvent<{
                eventData: any;
                args: OnBeforeEditCellEventArgs;
            }>,
        ) => {
            const { item, column } = e.detail.args;

            const currentOrderedColumns = orderedParamColumnsRef.current;
            const columnName = column.field;

            const editedCol = currentOrderedColumns.find((col) => col.field === columnName);

            // Not an editable column, return
            if (!editedCol || !editedCol.editable) {
                oldCellValueRef.current = null;
                return;
            }

            const oldValue = item[columnName];

            // Store the old value in the ref if it exists as 'String'
            if (oldValue !== null && oldValue !== undefined) {
                oldCellValueRef.current = String(oldValue);
            } else {
                oldCellValueRef.current = null;
            }
        },
        [orderedParamColumnsRef],
    );

    const handleCellChange = useCallback(
        (
            e: CustomEvent<{
                eventData: any;
                args: OnCellChangeEventArgs;
            }>,
        ) => {
            const { args } = e.detail;
            const { row, column, item } = args;

            const currentSortedIndices = sortedFilteredIndicesRef.current;
            const currentOrderedColumns = orderedParamColumnsRef.current;

            if (!currentSortedIndices || !currentOrderedColumns) {
                console.error("Values not loaded yet");
                setFeedbackAlert({
                    type: "warning",
                    message: "Values not loaded yet. Please wait...",
                    title: "Edit Warning",
                });
                return;
            }

            if (row < 0 || row >= currentSortedIndices.length) {
                console.error(`Row index ${row} is out of bounds`);
                setFeedbackAlert({
                    type: "error",
                    message: `Row index ${row} is out of bounds`,
                    title: "Edit Error",
                });
                return;
            }

            const columnName = column.field;
            const updatedValue = item[columnName];
            const dataIndex = currentSortedIndices[row];
            const editedCol = currentOrderedColumns.find((col) => col.field === columnName);
            const oldValue = oldCellValueRef.current;

            try {
                if (!editedCol) {
                    console.error("Column not found");
                    throw new Error("Column not found");
                }

                if (!editedCol.editable) {
                    console.error(`Column ${columnName} not editable`);
                    throw new Error(`Column ${columnName} not editable`);
                }

                const newValueString = updatedValue !== null && updatedValue !== undefined ? String(updatedValue) : "";

                // Set new value
                setCellValueFromString(editedCol, dataIndex, newValueString);

                setFeedbackAlert({
                    type: "success",
                    message: `Updated value ${oldValue} with ${updatedValue} in column: ${columnName}`,
                    title: "Edit Successful",
                });

                // Update the dataStore and rerender the grid
                dataStore.dataChanged([columnName]);
                const grid = gridRef.current?.slickGrid;
                if (grid) {
                    grid.invalidate();
                }
            } catch (err) {
                const error =
                    err instanceof Error ? err : new Error("An error occurred while trying to edit the value");
                setFeedbackAlert({
                    type: "error",
                    message: error.message,
                    title: "Edit Error",
                    stack: error.stack,
                    metadata: {
                        columnName,
                        oldValue,
                        newValue: updatedValue,
                    },
                });
            }
        },
        [dataStore, gridRef, orderedParamColumnsRef, sortedFilteredIndicesRef, setFeedbackAlert],
    );

    return {
        handleBeforeEditCell,
        handleCellChange,
    };
};

export default useEditCell;
