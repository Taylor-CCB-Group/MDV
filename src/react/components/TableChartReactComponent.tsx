import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import {
    type Column,
    Editors,
    type GridOption,
    type OnBeforeEditCellEventArgs,
    type OnCellChangeEventArgs,
    SlickgridReact,
    type SlickgridReactInstance,
} from "slickgrid-react";
import { useChartID, useConfig, useOrderedParamColumns } from "../hooks";
import type { TableChartReactConfig } from "./TableChartReactWrapper";
import SlickGridDataProvider from "../utils/SlickGridDataProvider";
import { useChart, useDataStore } from "../context";
import { runInAction } from "mobx";
import useSortedIndices from "../hooks/useSortedIndices";
import FindAndReplaceDialog from "./FindAndReplaceDialog";
import { useHighlightedIndex } from "../selectionHooks";
import type BaseChart from "@/charts/BaseChart";
import type { DataType, LoadedDataColumn } from "@/charts/charts";
import { replaceMatches } from "../utils/valueReplacementUtil";
import useFindReplace from "../hooks/useFindReplace";



// Type for chart with isFindReplaceOpen property
// type TableChartInstance = BaseChart<TableChartReactConfig> & {
//     isFindReplaceOpen: boolean;
// };

// function getColumnEditor(col: LoadedDataColumn<DataType>) {
//     if (!col.editable) return undefined; // not editable at all

//     switch (col.datatype) {
//         case "integer":
//         case "int32":
//             // return { model: Editors.integer };
//         case "double":
//             return { model: Editors.float };
//         default:
//             // basic text editor is fine, we’ll convert in our write logic
//             return { model: Editors.text };
//     }
// }

// todo: Add find and replace and editing functionality
// todo: Fix the styling of the table to have the scroll bars visible
const TableChartReactComponent = observer(() => {
    // Hooks
    const config = useConfig<TableChartReactConfig>();
    const dataStore = useDataStore();
    const chartId = useChartID();
    const chart = useChart<TableChartReactConfig>();
    const orderedParamColumns = useOrderedParamColumns();
    const sortedIndices = useSortedIndices();
    const highlightedIndex = useHighlightedIndex();

    // Find and replace states
    const [isFindReplaceOpen, setIsFindReplaceOpen] = useState(false);
    // const [foundMatches, setFoundMatches] = useState<FoundMatch[]>([]);
    const [searchColumn, setSearchColumn] = useState<string | null>(null);
    // const [matchCount, setMatchCount] = useState<number | null>(null);
    // const [currentMatchIndex, setCurrentMatchIndex] = useState(-1);

    // Refs
    const sortedIndicesRef = useRef(sortedIndices);
    const orderedParamColumnsRef = useRef(orderedParamColumns);
    const dataStoreRef = useRef(dataStore);
    const chartRef = useRef(chart);
    const oldCellValueRef = useRef<string | null>(null);
    const isSelectingRef = useRef(false); // Flag to prevent grid update during selection
    const gridRef = useRef<SlickgridReactInstance | null>(null);

    useEffect(() => {
        sortedIndicesRef.current = sortedIndices;
        orderedParamColumnsRef.current = orderedParamColumns;
        dataStoreRef.current = dataStore;
        chartRef.current = chart;
    }, [sortedIndices, dataStore, chart, orderedParamColumns]);

    const columnDefs = useMemo<Column[]>(() => {
        const cols: Column[] = [];

        if (config.include_index) {
            cols.push({
                id: "__index__",
                field: "__index__",
                name: "Index",
                sortable: true,
                minWidth: config.column_widths?.["__index__"] || 100,
            });
        }

        for (const col of orderedParamColumns) {
            // const editor = getColumnEditor(col);
            const isColumnEditable = col?.editable ?? false;
            cols.push({
                id: col.field,
                field: col.field,
                name: col.name,
                sortable: true,
                minWidth: config.column_widths?.[col.field] || 100,
                editor: isColumnEditable ? {model: Editors.text} : null,
                cssClass: isColumnEditable ? "mdv-editable-cell" : "",
                header: {
                    menu: {
                        commandItems: [
                            {
                                command: "find-replace",
                                title: "Find & Replace",
                                iconCssClass: "mdi mdi-magnify",
                            },
                        ],
                    },
                },
            });
        }

        return cols;
    }, [config.include_index, config.column_widths, orderedParamColumns]);

    const dataProvider = useMemo(() => {
        return new SlickGridDataProvider(dataStore, orderedParamColumns, sortedIndices, config.include_index);
    }, [dataStore, orderedParamColumns, sortedIndices, config.include_index]);

    const isColumnEditable = useMemo(() => {
        const column = orderedParamColumns.find((col) => col.field === searchColumn);
        return column?.editable ?? false;
    }, [searchColumn, orderedParamColumns]);

    useEffect(() => {
        const grid = gridRef.current?.slickGrid;
        if (grid && dataProvider) {
            // Skip update if we're in the middle of selection/navigation
            if (isSelectingRef.current) {
                console.log("Skipping grid update during selection");
                return;
            }
            grid.setData(dataProvider, true);
            grid.invalidate();
            console.log("Grid updated");
        }
    }, [dataProvider]);

    const handleGridCreated = useCallback(
        (e: CustomEvent<SlickgridReactInstance>) => {
            console.log("Grid created");
            gridRef.current = e.detail;

            const grid = e.detail.slickGrid;
            if (grid && dataProvider) {
                grid.setData(dataProvider, true);
                grid.render();
            }
        },
        [dataProvider],
    );

    useEffect(() => {
        const grid = gridRef.current?.slickGrid;
        if (!grid) return;

        const selectionHandler: any = grid.onSelectedRowsChanged.subscribe((_e, args) => {
            const selectedRows = args.rows;

            if (selectedRows.length > 0) {
                const indices = selectedRows.map((row) => sortedIndicesRef.current[row]);

                // Set flag to prevent grid update during selection
                isSelectingRef.current = true;

                dataStoreRef.current.dataHighlighted(indices, chartRef.current);

                // Reset flag after reactions settle
                setTimeout(() => {
                    isSelectingRef.current = false;
                }, 100);
            }
        });

        const sortHandler: any = grid.onSort.subscribe((_e, args) => {
            if ("sortCol" in args && args.sortCol && "sortAsc" in args) {
                const columnId = args.sortCol.field as string;
                const sortAsc = args.sortAsc as boolean;
                console.log("Sort event:", columnId, sortAsc ? "asc" : "desc");
                runInAction(() => {
                    config.sort = { columnId, ascending: sortAsc };
                });
            }
        });

        const headerMenuHandler = grid
            .getPubSubService()
            ?.subscribe("onHeaderMenuCommand", (event: { column: Column; command: string }) => {
                const { column, command } = event;
                // Remove Sort
                if (command === "clear-sort") {
                    console.log("Clear sort");
                    runInAction(() => {
                        config.sort = undefined;
                    });
                    // Sort Ascending
                } else if (command === "sort-asc") {
                    console.log("Sort Ascending");
                    runInAction(() => {
                        config.sort = { columnId: column.field, ascending: true };
                    });
                    // Sort Descending
                } else if (command === "sort-desc") {
                    console.log("Sort Descending");
                    runInAction(() => {
                        config.sort = { columnId: column.field, ascending: false };
                    });
                    // Find and replace
                } else if (command === "find-replace") {
                    console.log("Find and Replace");
                    setSearchColumn(column.field);
                    setIsFindReplaceOpen(true);
                }
            });

        const gridMenuHandler = grid
            .getPubSubService()
            ?.subscribe("onGridMenuCommand", ({ command }: { command: string }) => {
                if (command === "clear-sorting") {
                    console.log("Clear All sort");
                    runInAction(() => {
                        config.sort = undefined;
                    });
                }
            });

        return () => {
            selectionHandler?.unsubscribe();
            sortHandler?.unsubscribe();
            headerMenuHandler?.unsubscribe();
            gridMenuHandler?.unsubscribe();
        };
    }, [config]);

    // Handle highlighted data from other charts
    useEffect(() => {
        if (highlightedIndex === -1) return;

        const grid = gridRef.current?.slickGrid;
        if (!grid) return;

        // Find the row position in our filtered/sorted indices
        const pos = sortedIndicesRef.current.indexOf(highlightedIndex);

        if (pos !== -1) {
            // Set flag to prevent triggering our own highlight event
            isSelectingRef.current = true;

            grid.scrollRowIntoView(pos, false);
            grid.setSelectedRows([pos]);

            setTimeout(() => {
                isSelectingRef.current = false;
            }, 100);
        }
    }, [highlightedIndex]);

    const {
        matchCount,
        disableFindPrev,
        disableFindNext,
        handleFind,
        handleFindNext,
        handleFindPrev,
        handleReplace,
        handleReplaceAll,
        onReset,
    } = useFindReplace(orderedParamColumns,
        sortedIndices,
        dataStore,
        searchColumn,
        config,
        gridRef,
        isSelectingRef);


    // const handleFind = useCallback(
    //     (findText: string) => {
    //         console.log("Find: ", findText);

    //         if (!findText || findText.trim() === "") {
    //             setFoundMatches([]);
    //             setCurrentMatchIndex(-1);
    //             setMatchCount(null);
    //             return;
    //         }

    //         const matches: FoundMatch[] = [];

    //         const colIndex = orderedParamColumns.findIndex((col) => col.field === searchColumn);

    //         if (colIndex === -1) {
    //             console.log("Column index not found:", searchColumn);
    //             return;
    //         }

    //         const column = orderedParamColumns[colIndex];

    //         if (!column) {
    //             console.log("Column not found:", searchColumn);
    //             return;
    //         }

    //         for (let row = 0; row < sortedIndices.length; row++) {
    //             const dataIndex = sortedIndices[row];

    //             const value = dataStore.getRowText(dataIndex, column.field);

    //             if (value) {
    //                 //? Should this be case-insensitive and partial match?
    //                 if (typeof value === "string" && value.toLowerCase().includes(findText.toLowerCase())) {
    //                     matches.push({
    //                         rowIndex: row,
    //                         dataIndex,
    //                         column: column.field,
    //                         columnIndex: colIndex + (config.include_index ? 1 : 0),
    //                         value,
    //                     });
    //                 } else if (typeof value === "number" && value === Number(findText)) {
    //                     matches.push({
    //                         rowIndex: row,
    //                         dataIndex,
    //                         column: column.field,
    //                         columnIndex: colIndex + (config.include_index ? 1 : 0),
    //                         value,
    //                     });
    //                 }
    //             }
    //         }

    //         console.log(`Found ${matches.length} matches`);
    //         setFoundMatches(matches);
    //         setMatchCount(matches.length);

    //         if (matches.length > 0) {
    //             setCurrentMatchIndex(0);
    //             const grid = gridRef.current?.slickGrid;
    //             if (grid) {
    //                 const match = matches[0];
    //                 // Set flag before navigation
    //                 isSelectingRef.current = true;

    //                 // grid.scrollRowIntoView(match.rowIndex, false);
    //                 // grid.scrollCellIntoView(match.rowIndex, match.columnIndex, false);
    //                 // grid.setActiveCell(match.rowIndex, match.columnIndex);

    //                 grid.gotoCell(match.rowIndex, match.columnIndex, false);
    //                 // Reset flag after delay
    //                 setTimeout(() => {
    //                     isSelectingRef.current = false;
    //                 }, 100);
    //             }
    //         } else {
    //             setCurrentMatchIndex(-1);
    //         }
    //     },
    //     [sortedIndices, orderedParamColumns, dataStore, config.include_index, searchColumn],
    // );

    // const handleFindNext = useCallback(() => {
    //     // Sanity checks
    //     if (foundMatches.length === 0) return;
    //     if (currentMatchIndex + 1 >= foundMatches.length) return;

    //     const nextIndex = currentMatchIndex + 1;
    //     setCurrentMatchIndex(nextIndex);
    //     const grid = gridRef.current?.slickGrid;
    //     const match = foundMatches[nextIndex];

    //     if (grid && match) {
    //         // Set flag before navigation
    //         isSelectingRef.current = true;

    //         grid.gotoCell(match.rowIndex, match.columnIndex, false);

    //         // Reset flag after delay
    //         setTimeout(() => {
    //             isSelectingRef.current = false;
    //         }, 100);
    //     }

    //     console.log(`Match ${nextIndex + 1} of ${foundMatches.length}`);
    // }, [foundMatches, currentMatchIndex]);

    // const handleFindPrev = useCallback(() => {
    //     // Sanity checks
    //     if (foundMatches.length === 0) return;
    //     if (currentMatchIndex === 0) return;

    //     const prevIndex = (currentMatchIndex - 1) % foundMatches.length;
    //     setCurrentMatchIndex(prevIndex);

    //     const grid = gridRef.current?.slickGrid;
    //     const match = foundMatches[prevIndex];

    //     if (grid && match) {
    //         // Set flag before navigation
    //         isSelectingRef.current = true;

    //         grid.gotoCell(match.rowIndex, match.columnIndex, false);

    //         // Reset flag after delay
    //         setTimeout(() => {
    //             isSelectingRef.current = false;
    //         }, 100);
    //     }

    //     console.log(`Match ${prevIndex + 1} of ${foundMatches.length}`);
    // }, [foundMatches, currentMatchIndex]);

    // const getValueIndex = useCallback((replaceValue: string, values: string[], maxValues: number) => {
    //     let valueIndex = values.indexOf(replaceValue);

    //     if (valueIndex === -1) {
    //         if (values.length >= maxValues) {
    //             throw new Error(`Column exceeded ${maxValues} values while adding: ${replaceValue}`);
    //         }
    //         values.push(replaceValue);
    //         valueIndex = values.length - 1;
    //     }

    //     return valueIndex;
    // }, []);

    // const replaceMatches = useCallback(
    //     (
    //         column: LoadedDataColumn<DataType> | undefined,
    //         findValue: string,
    //         replaceValue: string,
    //         dataIndex: number,
    //     ) => {
    //         if (!column) {
    //             console.error("No column found for replace value: ", replaceValue);
    //             return false;
    //         }

    //         if (!column.data) {
    //             console.error("No data found in the column: ", searchColumn);
    //             return false;
    //         }

    //         if (column.datatype === "text" || column.datatype === "text16") {
    //             if (!column.values) {
    //                 console.error("No values found in the column: ", searchColumn);
    //                 return false;
    //             }
    //             // Based on datasource.md
    //             const maxValues = column.datatype === "text" ? 256 : 65536;
    //             const valueIndex = getValueIndex(replaceValue, column.values, maxValues);
    //             column.data[dataIndex] = valueIndex;
    //             return true;
    //         }

    //         if (column.datatype === "double" || column.datatype === "int32" || column.datatype === "integer") {
    //             const replaceNumber = Number.parseFloat(replaceValue);
    //             console.log("replace value: ", replaceValue, typeof replaceValue);
    //             console.log("replace number: ", replaceNumber, typeof replaceNumber);
    //             if (!Number.isFinite(replaceNumber)) {
    //                 console.error("Non-numeric value: ", searchColumn);
    //                 return false;
    //             }

    //             //? Can we assign a decimal to the datatype integer or int32?
    //             //? What if a user inputs a decimal for those columns?

    //             column.data[dataIndex] = replaceNumber;
    //             return true;
    //         }

    //         if (column.datatype === "multitext") {
    //             const { values, data, stringLength } = column;
    //             if (!values || !stringLength) {
    //                 console.error("No values found in the column: ", searchColumn);
    //                 return false;
    //             }

    //             const baseIndex = dataIndex * stringLength;

    //             const findLower = findValue.toLowerCase();
    //             const replaceIndex = getValueIndex(replaceValue, column.values, 65536);

    //             let replaced = false;
    //             for (let rowIndex = 0; rowIndex < stringLength; rowIndex++) {
    //                 const index = data[baseIndex + rowIndex];

    //                 if (index === 65535) continue;

    //                 const currentValue = values[index];
    //                 if (currentValue.toLowerCase() === findLower) {
    //                     data[baseIndex + rowIndex] = replaceIndex;
    //                     replaced = true;
    //                 }
    //             }

    //             return replaced;
    //         }

    //         if (column.datatype === "unique") {
    //             const { data, stringLength } = column;

    //             if (!stringLength) {
    //                 console.error("Missing string length for the column: ", searchColumn);
    //                 return false;
    //             }

    //             const encoder = new TextEncoder();
    //             const replaceEncoded = encoder.encode(replaceValue);
    //             console.log("encoded replace value: ", replaceEncoded);
    //             if (replaceEncoded.length > stringLength) {
    //                 console.error("Value too long for the column: ", searchColumn);
    //                 return false;
    //             }

    //             const baseIndex = dataIndex * stringLength;

    //             if (baseIndex + stringLength > data.length) {
    //                 console.error("Data index out of bounds for the column: ", searchColumn);
    //                 return false;
    //             }

    //             console.log("baseIndex: ", baseIndex);
    //             console.log("string length: ", stringLength);

    //             for (let i = 0; i < stringLength; i++) {
    //                 if (i < replaceEncoded.length) {
    //                     console.log("in if value: ", data[baseIndex + i], replaceEncoded[i]);
    //                     data[baseIndex + i] = replaceEncoded[i];
    //                 } else {
    //                     console.log("in else value: ", data[baseIndex + i]);
    //                     data[baseIndex + i] = 0;
    //                 }
    //             }

    //             console.log("updated data: ", data);

    //             return true;
    //         }

    //         console.error("Replace not supported for datatype: ", column.datatype);
    //         return false;
    //     },
    //     [searchColumn, getValueIndex],
    // );

    // const handleReplace = useCallback(
    //     (findValue: string, replaceValue: string) => {
    //         if (!searchColumn) {
    //             console.log("No column selected for replace");
    //             return;
    //         }

    //         const column = orderedParamColumns.find((col) => col.field === searchColumn);

    //         if (!column?.editable) {
    //             console.error(`Column ${searchColumn} is not editable`);
    //             return;
    //         }

    //         if (foundMatches.length === 0 || currentMatchIndex < 0 || currentMatchIndex >= foundMatches.length) {
    //             console.log("No current match found for replace value: ", replaceValue);
    //             return;
    //         }

    //         const match = foundMatches[currentMatchIndex];

    //         const matchReplaced = replaceMatches(searchColumn, column, findValue, replaceValue, match.dataIndex);

    //         if (matchReplaced) {
    //             dataStore.dataChanged([searchColumn]);
    //             const grid = gridRef.current?.slickGrid;
    //             if (grid) {
    //                 grid.invalidate();
    //                 grid.render();
    //             }
    //             handleFind(findValue);
    //         }
    //     },
    //     [searchColumn, foundMatches, orderedParamColumns, currentMatchIndex, dataStore, handleFind],
    // );

    // const handleReplaceAll = useCallback(
    //     (findValue: string, replaceValue: string) => {
    //         if (!searchColumn) {
    //             console.log("No column selected for replace all");
    //             return;
    //         }

    //         const column = orderedParamColumns.find((col) => col.field === searchColumn);

    //         if (!column?.editable) {
    //             console.error(`Column ${searchColumn} is not editable`);
    //             return;
    //         }

    //         if (foundMatches.length === 0) {
    //             console.log("No matches to replace");
    //             return;
    //         }

    //         let valuesReplaced = false;
    //         for (const match of foundMatches) {
    //             const matchReplaced = replaceMatches(searchColumn, column, findValue, replaceValue, match.dataIndex);
    //             valuesReplaced = valuesReplaced || matchReplaced;
    //         }

    //         if (valuesReplaced) {
    //             dataStore.dataChanged([searchColumn]);
    //             const grid = gridRef.current?.slickGrid;
    //             if (grid) {
    //                 grid.invalidate();
    //                 grid.render();
    //             }
    //             handleFind(findValue);
    //         }
    //     },
    //     [searchColumn, foundMatches, orderedParamColumns, dataStore, handleFind],
    // );

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

            if (oldValue !== null || oldValue !== undefined) {
                console.log("old value: ", typeof oldValue);
                oldCellValueRef.current = String(oldValue);
            }
        },
        [],
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
        [dataStore],
    );

    const onClose = useCallback(() => {
        setIsFindReplaceOpen(false);
        setSearchColumn(null);
        // setFoundMatches([]);
        // setMatchCount(null);
        // setCurrentMatchIndex(-1);
        onReset();
    }, [onReset]);

    const options: GridOption = useMemo(
        () => ({
            gridWidth: "100%",
            gridHeight: "600px",
            darkMode: window?.mdv?.chartManager?.theme === "dark",
            enableSorting: true,
            multiColumnSort: false,
            enableGridMenu: false,
            enableHeaderMenu: true,
            alwaysShowVerticalScroll: true,
            alwaysAllowHorizontalScroll: true,
            enableAutoResize: true,
            autoResize: {
                container: `#react-grid-${chartId}`,
                calculateAvailableSizeBy: "container",
                resizeDetection: "container",
                autoHeight: true,
            },
            editable: true,
            enableCellNavigation: true,
            enableExcelCopyBuffer: true,
            multiSelect: true,
            enableRowSelection: true,
            rowSelectionOptions: {
                selectActiveRow: true,
            },
        }),
        [chartId],
    );

    return (
        <div id={`react-grid-${chartId}`} className="absolute w-[100%] h-[100%] slickgrid-react-container">
            <SlickgridReact
                gridId={`table-${chartId}`}
                columns={columnDefs}
                dataset={[]}
                options={options}
                onReactGridCreated={handleGridCreated}
                onBeforeEditCell={handleBeforeEditCell}
                onCellChange={handleCellChange}
            />

            <FindAndReplaceDialog
                open={isFindReplaceOpen}
                onClose={onClose}
                handleFind={handleFind}
                handleFindPrev={handleFindPrev}
                handleFindNext={handleFindNext}
                handleReplace={handleReplace}
                handleReplaceAll={handleReplaceAll}
                foundMatches={matchCount}
                columnName={searchColumn}
                // todo: might take a look at this to clean up
                // disableFindPrev={foundMatches.length === 0 || (foundMatches.length > 0 && currentMatchIndex <= 0)}
                // disableFindNext={
                //     foundMatches.length === 0 ||
                //     (foundMatches.length > 0 && currentMatchIndex + 1 >= foundMatches.length)
                // }
                disableFindPrev={disableFindPrev}
                disableFindNext={disableFindNext}
                isColumnEditable={isColumnEditable}
            />
        </div>
    );
});

export default TableChartReactComponent;
