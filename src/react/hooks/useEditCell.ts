import type { DataType, LoadedDataColumn } from "@/charts/charts";
import { useCallback, useRef } from "react";
import type { OnBeforeEditCellEventArgs, OnCellChangeEventArgs, SlickgridReactInstance } from "slickgrid-react";
import { replaceMatches } from "../utils/valueReplacementUtil";

const useEditCell = (
    orderedParamColumnsRef: React.MutableRefObject<LoadedDataColumn<DataType>[]>,
    sortedIndicesRef: React.MutableRefObject<Uint32Array>,
    dataStore: any,
    gridRef: React.MutableRefObject<SlickgridReactInstance | null>,
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
                console.log("old value: ", typeof oldValue);
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
                console.error("Values not available yet");
                return;
            }

            const columnName = column.field;
            const updatedValue = item[columnName];
            const dataIndex = currentSortedIndices[row];
            const editedCol = currentOrderedColumns.find((col) => col.field === columnName);
            const oldValue = oldCellValueRef.current;

            if (!editedCol || !editedCol?.editable) {
                return;
            }
            console.log("params: ", editedCol, oldValue, updatedValue, dataIndex, typeof oldValue, typeof updatedValue);
            const changed = replaceMatches(columnName, editedCol, oldValue as string, updatedValue, dataIndex);
            console.log("changed value", changed);
            if (changed) {
                dataStore.dataChanged([columnName]);
                const grid = gridRef.current?.slickGrid;
                if (grid) {
                    grid.invalidate();
                    grid.render();
                }
            }

        },
        [dataStore, gridRef, orderedParamColumnsRef, sortedIndicesRef],
    );

    return {
        handleBeforeEditCell,
        handleCellChange,
    };
};

export default useEditCell;