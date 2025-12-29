import type { DataType, LoadedDataColumn } from "@/charts/charts";
import { useCallback, useRef } from "react";
import type { OnBeforeEditCellEventArgs, OnCellChangeEventArgs, SlickgridReactInstance } from "slickgrid-react";
import { replaceMatches, setCellValueFromString } from "../utils/valueReplacementUtil";
import type DataStore from "@/datastore/DataStore";
import type { FeedbackAlert } from "../components/TableChartReactComponent";

/**
 *
 * Custom hook to handle editing of cell in the grid
 */
const useEditCell = (
    orderedParamColumnsRef: React.MutableRefObject<LoadedDataColumn<DataType>[]>,
    sortedIndicesRef: React.MutableRefObject<Uint32Array>,
    dataStore: DataStore,
    gridRef: React.MutableRefObject<SlickgridReactInstance | null>,
    setFeedbackAlert: (alert: FeedbackAlert) => void,
) => {
    const oldCellValueRef = useRef<string | null>(null);

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

            if (!editedCol || !editedCol.editable) {
                oldCellValueRef.current = null;
                return;
            }

            const oldValue = item[columnName];

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

            const currentSortedIndices = sortedIndicesRef.current;
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

                if (!editedCol?.editable) {
                    console.error(`Column ${columnName} not editable`);
                    throw new Error(`Column ${columnName} not editable`);
                }

                const newValueString = updatedValue !== null && updatedValue !== undefined ? String(updatedValue) : "";

                // setCellValueFromString now throws on error instead of returning false
                setCellValueFromString(editedCol, dataIndex, newValueString);

                // If we get here, it succeeded
                setFeedbackAlert({
                    type: "success",
                    message: `Updated value ${oldValue} with ${updatedValue} in column: ${columnName}`,
                    title: "Edit Successful",
                });
                dataStore.dataChanged([columnName]);
                const grid = gridRef.current?.slickGrid;
                if (grid) {
                    grid.invalidate();
                    grid.render();
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
        [dataStore, gridRef, orderedParamColumnsRef, sortedIndicesRef, setFeedbackAlert],
    );

    return {
        handleBeforeEditCell,
        handleCellChange,
    };
};

export default useEditCell;
