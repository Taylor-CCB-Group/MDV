import type { DataType, LoadedDataColumn } from "@/charts/charts";
import { useCallback, useMemo, useState } from "react";
import { replaceMatches } from "../utils/valueReplacementUtil";
import type { SlickgridReactInstance } from "slickgrid-react";

export type FoundMatch = {
    value: string | number;
    column: string;
    columnIndex: number;
    rowIndex: number;
    dataIndex: number;
};

const useFindReplace = (
    orderedParamColumns: LoadedDataColumn<DataType>[],
    sortedIndices: Uint32Array,
    dataStore: any,
    searchColumn: string | null,
    config: { include_index?: boolean },
    gridRef: React.MutableRefObject<SlickgridReactInstance | null>,
    isSelectingRef: React.MutableRefObject<boolean>,
) => {
    const [foundMatches, setFoundMatches] = useState<FoundMatch[]>([]);
    const [matchCount, setMatchCount] = useState<number | null>(null);
    const [currentMatchIndex, setCurrentMatchIndex] = useState(-1);

    const disableFindPrev = useMemo(() => {
        return foundMatches.length === 0 || (foundMatches.length > 0 && currentMatchIndex <= 0);
    }, [currentMatchIndex, foundMatches.length]);

    const disableFindNext = useMemo(() => {
        return foundMatches.length === 0 || (foundMatches.length > 0 && currentMatchIndex + 1 >= foundMatches.length);
    }, [currentMatchIndex, foundMatches.length]);

    const handleFind = useCallback(
        (findText: string) => {
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
                return;
            }

            const column = orderedParamColumns[colIndex];

            if (!column) {
                console.log("Column not found:", searchColumn);
                return;
            }

            for (let row = 0; row < sortedIndices.length; row++) {
                const dataIndex = sortedIndices[row];

                const value = dataStore.getRowText(dataIndex, column.field);

                if (value) {
                    //? Should this be case-insensitive and partial match?
                    if (typeof value === "string" && value.toLowerCase().includes(findText.toLowerCase())) {
                        matches.push({
                            rowIndex: row,
                            dataIndex,
                            column: column.field,
                            columnIndex: colIndex + (config.include_index ? 1 : 0),
                            value,
                        });
                    } else if (typeof value === "number" && value === Number(findText)) {
                        matches.push({
                            rowIndex: row,
                            dataIndex,
                            column: column.field,
                            columnIndex: colIndex + (config.include_index ? 1 : 0),
                            value,
                        });
                    }
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
                    // Reset flag after delay
                    setTimeout(() => {
                        isSelectingRef.current = false;
                    }, 100);
                }
            } else {
                setCurrentMatchIndex(-1);
            }
        },
        [sortedIndices, orderedParamColumns, dataStore, config.include_index, searchColumn, gridRef, isSelectingRef],
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

            // Reset flag after delay
            setTimeout(() => {
                isSelectingRef.current = false;
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

            // Reset flag after delay
            setTimeout(() => {
                isSelectingRef.current = false;
            }, 100);
        }

        console.log(`Match ${prevIndex + 1} of ${foundMatches.length}`);
    }, [foundMatches, currentMatchIndex, gridRef, isSelectingRef]);

    const handleReplace = useCallback(
        (findValue: string, replaceValue: string) => {
            if (!searchColumn) {
                console.log("No column selected for replace");
                return;
            }

            const column = orderedParamColumns.find((col) => col.field === searchColumn);

            if (!column?.editable) {
                console.error(`Column ${searchColumn} is not editable`);
                return;
            }

            if (foundMatches.length === 0 || currentMatchIndex < 0 || currentMatchIndex >= foundMatches.length) {
                console.log("No current match found for replace value: ", replaceValue);
                return;
            }

            const match = foundMatches[currentMatchIndex];

            const matchReplaced = replaceMatches(searchColumn, column, findValue, replaceValue, match.dataIndex);

            if (matchReplaced) {
                dataStore.dataChanged([searchColumn]);
                const grid = gridRef.current?.slickGrid;
                if (grid) {
                    grid.invalidate();
                    grid.render();
                }
                handleFind(findValue);
            }
        },
        [searchColumn, foundMatches, orderedParamColumns, currentMatchIndex, dataStore, handleFind, gridRef],
    );

    const handleReplaceAll = useCallback(
        (findValue: string, replaceValue: string) => {
            if (!searchColumn) {
                console.log("No column selected for replace all");
                return;
            }

            const column = orderedParamColumns.find((col) => col.field === searchColumn);

            if (!column?.editable) {
                console.error(`Column ${searchColumn} is not editable`);
                return;
            }

            if (foundMatches.length === 0) {
                console.log("No matches to replace");
                return;
            }

            let valuesReplaced = false;
            for (const match of foundMatches) {
                const matchReplaced = replaceMatches(searchColumn, column, findValue, replaceValue, match.dataIndex);
                valuesReplaced = valuesReplaced || matchReplaced;
            }

            if (valuesReplaced) {
                dataStore.dataChanged([searchColumn]);
                const grid = gridRef.current?.slickGrid;
                if (grid) {
                    grid.invalidate();
                    grid.render();
                }
                handleFind(findValue);
            }
        },
        [searchColumn, foundMatches, orderedParamColumns, dataStore, handleFind, gridRef],
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
