import type { DataType, LoadedDataColumn } from "@/charts/charts";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { replaceMatches, replaceValueInString } from "../utils/valueReplacementUtil";
import type { SlickgridReactInstance } from "slickgrid-react";
import type DataStore from "@/datastore/DataStore";
import type { FeedbackAlert } from "../components/TableChartReactComponent";

export type FoundMatch = {
    value: string | number;
    column: string;
    columnIndex: number;
    rowIndex: number;
    dataIndex: number;
};

/**
 * Custom hook which handles find and replace logic
 */
const useFindReplace = (
    orderedParamColumns: LoadedDataColumn<DataType>[],
    sortedFilteredIndices: Uint32Array,
    dataStore: DataStore,
    searchColumn: string | null,
    config: { include_index?: boolean },
    gridRef: React.MutableRefObject<SlickgridReactInstance | null>,
    isSelectingRef: React.MutableRefObject<boolean>,
    setFeedbackAlert: (alert: FeedbackAlert) => void,
) => {
    const [foundMatches, setFoundMatches] = useState<FoundMatch[]>([]);
    const [matchCount, setMatchCount] = useState<number | null>(null);
    const [currentMatchIndex, setCurrentMatchIndex] = useState(-1);

    const selectionTimeoutRef = useRef<ReturnType<typeof setTimeout> | null>(null);


    useEffect(() => {
        return () => {
            if (selectionTimeoutRef.current !== null) {
                clearTimeout(selectionTimeoutRef.current);
            }
        };
    }, []);

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

                if (matches.length > 0) {
                    setCurrentMatchIndex(0);
                    const grid = gridRef.current?.slickGrid;
                    if (grid) {
                        const match = matches[0];
                        // Set flag before navigation
                        isSelectingRef.current = true;

                        grid.gotoCell(match.rowIndex, match.columnIndex, false);

                        // Clear any existing timeout before scheduling a new one
                        if (selectionTimeoutRef.current !== null) {
                            clearTimeout(selectionTimeoutRef.current);
                        }

                        // Reset flag after delay
                        selectionTimeoutRef.current = setTimeout(() => {
                            isSelectingRef.current = false;
                            selectionTimeoutRef.current = null;
                        }, 100);
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
            isSelectingRef,
            setFeedbackAlert,
        ],
    );

    const handleFindNext = useCallback(() => {
        // Sanity checks
        if (foundMatches.length === 0) return;
        if (currentMatchIndex + 1 >= foundMatches.length) return;

        const nextIndex = currentMatchIndex + 1;
        setCurrentMatchIndex(nextIndex);
        const grid = gridRef.current?.slickGrid;
        const match = foundMatches[nextIndex];

        if (grid && match) {
            // Set flag before navigation
            isSelectingRef.current = true;

            grid.gotoCell(match.rowIndex, match.columnIndex, false);

            // Clear any existing timeout before scheduling a new one
            if (selectionTimeoutRef.current !== null) {
                clearTimeout(selectionTimeoutRef.current);
            }

            // Reset flag after delay
            selectionTimeoutRef.current = setTimeout(() => {
                isSelectingRef.current = false;
                selectionTimeoutRef.current = null;
            }, 100);
        }

        console.log(`Match ${nextIndex + 1} of ${foundMatches.length}`);
    }, [foundMatches, currentMatchIndex, gridRef, isSelectingRef]);

    const handleFindPrev = useCallback(() => {
        // Sanity checks
        if (foundMatches.length === 0) return;
        if (currentMatchIndex === 0) return;

        const prevIndex = (currentMatchIndex - 1) % foundMatches.length;
        setCurrentMatchIndex(prevIndex);

        const grid = gridRef.current?.slickGrid;
        const match = foundMatches[prevIndex];

        if (grid && match) {
            // Set flag before navigation
            isSelectingRef.current = true;

            grid.gotoCell(match.rowIndex, match.columnIndex, false);

            // Clear any existing timeout before scheduling a new one
            if (selectionTimeoutRef.current !== null) {
                clearTimeout(selectionTimeoutRef.current);
            }

            // Reset flag after delay
            selectionTimeoutRef.current = setTimeout(() => {
                isSelectingRef.current = false;
                selectionTimeoutRef.current = null;
            }, 100);
        }

        console.log(`Match ${prevIndex + 1} of ${foundMatches.length}`);
    }, [foundMatches, currentMatchIndex, gridRef, isSelectingRef]);

    //! Need a check for replace value
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

                const matchReplaced = replaceMatches(searchColumn, column, findValue, replaceValue, match.dataIndex);

                if (matchReplaced) {
                    setFeedbackAlert({
                        type: "success",
                        message: `Replaced ${findValue} with ${replaceValue} in column: ${searchColumn}`,
                        title: "Replace Successful",
                    });
                    dataStore.dataChanged([searchColumn]);
                    const grid = gridRef.current?.slickGrid;
                    if (grid) {
                        grid.invalidate();
                        grid.render();
                    }
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

                // let valuesReplaced = false;
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

                if (successCount > 0) {
                    dataStore.dataChanged([searchColumn]);
                    const grid = gridRef.current?.slickGrid;
                    if (grid) {
                        grid.invalidate();
                        grid.render();
                    }
                    handleFind(findValue);
                }

                if (successCount > 0 && !errorCount) {
                    setFeedbackAlert({
                        type: "success",
                        message: `Replaced ${successCount} occurrences of ${findValue} with ${replaceValue} in column: ${searchColumn}`,
                        title: "Replace All Successful",
                    });
                    return;
                } else if (!successCount && !errorCount) {
                    setFeedbackAlert({
                        type: "warning",
                        message: `No matches were replaced for ${findValue} with ${replaceValue} in column: ${searchColumn}`,
                        title: "Replace All Warning",
                    });
                    return;
                } else {
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
