import type { DataType, LoadedDataColumn } from "@/charts/charts";
import { useCallback, useMemo, useState } from "react";
import { replaceMatches } from "../utils/valueReplacementUtil";
import type { SlickgridReactInstance } from "slickgrid-react";
import type DataStore from "@/datastore/DataStore";
import type { FeedbackAlert } from "../components/FeedbackAlertComponent";

export type FoundMatch = {
    value: string | number;
    column: string;
    columnIndex: number;
    rowIndex: number;
    dataIndex: number;
};

/**
 * Hook for handling find and replace functionality
 * 
 * Search runs over the currently selected column and uses sortedFilteredIndices (visible rows).
 * Replace uses valueReplacementUtil and notifies the data store so changes persist.
 */
const useFindReplace = (
    orderedParamColumns: LoadedDataColumn<DataType>[],
    sortedFilteredIndices: Uint32Array,
    dataStore: DataStore,
    searchColumn: string | null,
    config: { include_index?: boolean },
    gridRef: React.MutableRefObject<SlickgridReactInstance | null>,
    selectionSourceRef: React.MutableRefObject<'user' | 'programmatic' | null>,
    setFeedbackAlert: (alert: FeedbackAlert) => void,
) => {
    const [foundMatches, setFoundMatches] = useState<FoundMatch[]>([]);
    const [matchCount, setMatchCount] = useState<number | null>(null);
    const [currentMatchIndex, setCurrentMatchIndex] = useState(-1);

    // Disable flags for the find and replace dialog to disable the prev and next buttons
    const disableFindPrev = useMemo(() => {
        return foundMatches.length === 0 || (foundMatches.length > 0 && currentMatchIndex <= 0);
    }, [currentMatchIndex, foundMatches.length]);

    const disableFindNext = useMemo(() => {
        return foundMatches.length === 0 || (foundMatches.length > 0 && currentMatchIndex + 1 >= foundMatches.length);
    }, [currentMatchIndex, foundMatches.length]);

    const handleFind = useCallback(
        (findText: string) => {
            try {
                console.log("Find: ", findText);

                if (!findText || findText.trim() === "") {
                    setFoundMatches([]);
                    setCurrentMatchIndex(-1);
                    setMatchCount(null);
                    return;
                }

                const matches: FoundMatch[] = [];

                const colIndex = orderedParamColumns.findIndex((col) => col.field === searchColumn);

                if (colIndex === -1) {
                    console.log("Column index not found:", searchColumn);
                    throw new Error(`Column index not found: ${searchColumn}`);
                }

                const column = orderedParamColumns[colIndex];

                if (!column) {
                    console.log("Column not found:", searchColumn);
                    throw new Error(`Column not found: ${searchColumn}`);
                }

                // Find the matches in the sortedFilteredIndices array
                for (let row = 0; row < sortedFilteredIndices.length; row++) {
                    const dataIndex = sortedFilteredIndices[row];

                    const value = column.getValue(dataIndex);

                    if (value === null) continue;

                    const displayedValue = String(value);
                    const lowerDisplayedValue = displayedValue.toLowerCase();
                    
                    if (lowerDisplayedValue.includes(findText.toLowerCase())) {
                        matches.push({
                            rowIndex: row,
                            dataIndex,
                            column: column.field,
                            columnIndex: colIndex + (config.include_index ? 1 : 0),
                            value: displayedValue,
                        });
                    }
                }

                console.log(`Found ${matches.length} matches`);
                setFoundMatches(matches);
                setMatchCount(matches.length);

                // If matches are found, set match index to 0 and navigate to the cell
                if (matches.length > 0) {
                    setCurrentMatchIndex(0);
                    const grid = gridRef.current?.slickGrid;
                    if (grid) {
                        const match = matches[0];
                        selectionSourceRef.current = 'programmatic';
                        grid.gotoCell(match.rowIndex, match.columnIndex, false);
                        requestAnimationFrame(() => {
                            selectionSourceRef.current = null;
                        });
                    }
                } else {
                    setCurrentMatchIndex(-1);
                }
            } catch (err) {
                const error =
                    err instanceof Error
                        ? err
                        : new Error(`An error occurred while trying to find the value: ${findText}`);
                setFeedbackAlert({
                    type: "error",
                    message: error.message,
                    stack: error.stack,
                    metadata: {
                        columnName: searchColumn,
                        findValue: findText,
                    },
                    title: "Find Error",
                });
            }
        },
        [
            sortedFilteredIndices,
            orderedParamColumns,
            config.include_index,
            searchColumn,
            gridRef,
            selectionSourceRef,
            setFeedbackAlert,
        ],
    );

    const handleFindNext = useCallback(() => {
        // Sanity checks
        if (foundMatches.length === 0) return;
        if (currentMatchIndex + 1 >= foundMatches.length) return;

        // Update match index
        const nextIndex = currentMatchIndex + 1;
        setCurrentMatchIndex(nextIndex);
        const grid = gridRef.current?.slickGrid;
        const match = foundMatches[nextIndex];

        // Go to the cell in grid's display
        if (grid && match) {
            selectionSourceRef.current = 'programmatic';
            grid.gotoCell(match.rowIndex, match.columnIndex, false);
            requestAnimationFrame(() => {
                selectionSourceRef.current = null;
            });
        }

        console.log(`Match ${nextIndex + 1} of ${foundMatches.length}`);
    }, [foundMatches, currentMatchIndex, gridRef, selectionSourceRef]);

    const handleFindPrev = useCallback(() => {
        // Sanity checks
        if (foundMatches.length === 0) return;
        if (currentMatchIndex === 0) return;

        // Update match index
        const prevIndex = (currentMatchIndex - 1) % foundMatches.length;
        setCurrentMatchIndex(prevIndex);

        const grid = gridRef.current?.slickGrid;
        const match = foundMatches[prevIndex];

        // Go to the cell in grid's display
        if (grid && match) {
            selectionSourceRef.current = 'programmatic';
            grid.gotoCell(match.rowIndex, match.columnIndex, false);
            requestAnimationFrame(() => {
                selectionSourceRef.current = null;
            });
        }

        console.log(`Match ${prevIndex + 1} of ${foundMatches.length}`);
    }, [foundMatches, currentMatchIndex, gridRef, selectionSourceRef]);

    const handleReplace = useCallback(
        (findValue: string, replaceValue: string) => {
            try {
                if (!searchColumn) {
                    console.log("No column selected for replace");
                    throw new Error("No column selected for replace");
                }

                const column = orderedParamColumns.find((col) => col.field === searchColumn);

                if (!column?.editable) {
                    console.error(`Column ${searchColumn} is not editable`);
                    throw new Error(`Column ${searchColumn} is not editable`);
                }

                if (foundMatches.length === 0 || currentMatchIndex < 0 || currentMatchIndex >= foundMatches.length) {
                    console.log("No current match found for replace value: ", replaceValue);
                    throw new Error(`No current match found for replace value: ${replaceValue}`);
                }

                const match = foundMatches[currentMatchIndex];

                // Replace the cell value with new value
                const matchReplaced = replaceMatches(searchColumn, column, findValue, replaceValue, match.dataIndex);

                if (matchReplaced) {
                    setFeedbackAlert({
                        type: "success",
                        message: `Replaced ${findValue} with ${replaceValue} in column: ${searchColumn}`,
                        title: "Replace Successful",
                    });

                    // Update the data store and rerender the grid
                    dataStore.dataChanged([searchColumn]);
                    const grid = gridRef.current?.slickGrid;
                    if (grid) {
                        grid.invalidate();
                    }

                    // Find the rest of the occurrences
                    handleFind(findValue);
                } else {
                    // No match found - this is not an error, just inform the user
                    setFeedbackAlert({
                        type: "warning",
                        message: `No match found for "${findValue}" in the current cell`,
                        title: "Replace Warning",
                    });
                }
            } catch (err) {
                const error =
                    err instanceof Error ? err : new Error("An error occurred while trying to replace the value");
                setFeedbackAlert({
                    type: "error",
                    message: error.message,
                    stack: error.stack,
                    metadata: {
                        columnName: searchColumn,
                        findValue: findValue,
                        replaceValue: replaceValue,
                    },
                    title: "Replace Error",
                });
            }
        },
        [
            searchColumn,
            foundMatches,
            orderedParamColumns,
            currentMatchIndex,
            dataStore,
            handleFind,
            gridRef,
            setFeedbackAlert,
        ],
    );

    const handleReplaceAll = useCallback(
        (findValue: string, replaceValue: string) => {
            try {
                if (!searchColumn) {
                    console.log("No column selected for replace all");
                    throw new Error("No column selected for replace all");
                }

                const column = orderedParamColumns.find((col) => col.field === searchColumn);

                if (!column?.editable) {
                    console.error(`Column ${searchColumn} is not editable`);
                    throw new Error(`Column ${searchColumn} is not editable`);
                }

                if (foundMatches.length === 0) {
                    console.log("No match found for replace value: ", replaceValue);
                    throw new Error(`No match found for replace value: ${replaceValue}`);
                }

                let successCount = 0;
                let errorCount = 0;

                // Replace all matches with new value
                for (const match of foundMatches) {
                    try {
                        const matchReplaced = replaceMatches(
                            searchColumn,
                            column,
                            findValue,
                            replaceValue,
                            match.dataIndex,
                        );

                        if (matchReplaced) {
                            successCount++;
                        }
                    } catch (err) {
                        errorCount++;
                    }
                }

                // Update the data store and rerender the grid if values were changed
                if (successCount > 0) {
                    dataStore.dataChanged([searchColumn]);
                    const grid = gridRef.current?.slickGrid;
                    if (grid) {
                        grid.invalidate();
                    }
                    handleFind(findValue);
                }

                if (successCount > 0 && !errorCount) {
                    // Display success if there are no errors
                    setFeedbackAlert({
                        type: "success",
                        message: `Replaced ${successCount} occurrences of ${findValue} with ${replaceValue} in column: ${searchColumn}`,
                        title: "Replace All Successful",
                    });
                    return;
                } else if (!successCount && !errorCount) {
                    // Display warning if there was no replacement made
                    setFeedbackAlert({
                        type: "warning",
                        message: `No matches were replaced for ${findValue} with ${replaceValue} in column: ${searchColumn}`,
                        title: "Replace All Warning",
                    });
                    return;
                } else {
                    // Display error if there were one or more errors
                    setFeedbackAlert({
                        type: "error",
                        message: `${errorCount} error(s) occurred while trying to replace all occurrences of ${findValue} with ${replaceValue} in column: ${searchColumn}`,
                        title: "Replace All Error",
                        metadata: {
                            columnName: searchColumn,
                            findValue,
                            replaceValue,
                            errorCount,
                        },
                    });
                    return;
                }
            } catch (err) {
                const error =
                    err instanceof Error
                        ? err
                        : new Error(
                              `An error occurred while trying to replace all occurrences of ${findValue} with ${replaceValue}`,
                          );
                setFeedbackAlert({
                    type: "error",
                    message: error.message,
                    stack: error.stack,
                    metadata: {
                        columnName: searchColumn,
                        findValue: findValue,
                        replaceValue: replaceValue,
                    },
                    title: "Replace All Error",
                });
            }
        },
        [searchColumn, foundMatches, orderedParamColumns, dataStore, handleFind, gridRef, setFeedbackAlert],
    );

    const onReset = useCallback(() => {
        setFoundMatches([]);
        setMatchCount(null);
        setCurrentMatchIndex(-1);
    }, []);

    return {
        foundMatches,
        matchCount,
        currentMatchIndex,
        disableFindNext,
        disableFindPrev,
        setFoundMatches,
        setMatchCount,
        setCurrentMatchIndex,
        handleFind,
        handleFindNext,
        handleFindPrev,
        handleReplace,
        handleReplaceAll,
        onReset,
    };
};

export default useFindReplace;
