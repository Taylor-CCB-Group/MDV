import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { type Column, type GridOption, SlickgridReact, type SlickgridReactInstance } from "slickgrid-react";
import { useChartID, useConfig, useOrderedParamColumns } from "../hooks";
import type { TableChartReactConfig } from "./TableChartReactWrapper";
import SlickGridDataProvider from "../utils/SlickGridDataProvider";
import { useChart, useDataStore } from "../context";
import { runInAction } from "mobx";
import useSortedIndices from "../hooks/useSortedIndices";
import FindAndReplaceDialog from "./FindAndReplaceDialog";
import { useHighlightedIndex } from "../selectionHooks";
import type BaseChart from "@/charts/BaseChart";

export type FindMatch = {
    value: string | number;
    column: string;
    columnIndex: number;
    rowIndex: number;
    dataIndex: number;
};

// Type for chart with isFindReplaceOpen property
// type TableChartInstance = BaseChart<TableChartReactConfig> & {
//     isFindReplaceOpen: boolean;
// };

// todo: Add find and replace and editing functionality
// todo: Fix the styling of the table to have the scroll bars visible
const TableChartReactComponent = observer(() => {
    const config = useConfig<TableChartReactConfig>();
    const dataStore = useDataStore();
    const chartId = useChartID();
    const chart = useChart<TableChartReactConfig>();
    const orderedParamColumns = useOrderedParamColumns();
    const filteredIndices = useSortedIndices();
    const highlightedIndex = useHighlightedIndex();
    const [isFindReplaceOpen, setIsFindReplaceOpen] = useState(false);
    const [findMatches, setFindMatches] = useState<FindMatch[]>([]);
    const [matchCount, setMatchCount] = useState<number | null>(null);
    const [currentMatchIndex, setCurrentMatchIndex] = useState(-1);
    const filteredIndicesRef = useRef(filteredIndices);
    const dataStoreRef = useRef(dataStore);
    const chartRef = useRef(chart);
    const isSelectingRef = useRef(false); // Flag to prevent grid update during selection
    const [findColumn, setFindColumn] = useState<string | null>(null);
    const gridRef = useRef<SlickgridReactInstance | null>(null);

    useEffect(() => {
        filteredIndicesRef.current = filteredIndices;
        dataStoreRef.current = dataStore;
        chartRef.current = chart;
    }, [filteredIndices, dataStore, chart]);

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
            cols.push({
                id: col.field,
                field: col.field,
                name: col.name,
                sortable: true,
                minWidth: config.column_widths?.[col.field] || 100,
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
        return new SlickGridDataProvider(dataStore, orderedParamColumns, filteredIndices, config.include_index);
    }, [dataStore, orderedParamColumns, filteredIndices, config.include_index]);

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
                const indices = selectedRows.map((row) => filteredIndicesRef.current[row]);

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
                    setFindColumn(column.field);
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
        const pos = filteredIndicesRef.current.indexOf(highlightedIndex);

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

    const handleFind = useCallback(
        (findText: string) => {
            console.log("Find: ", findText);

            if (!findText || findText.trim() === "") {
                setFindMatches([]);
                setCurrentMatchIndex(-1);
                setMatchCount(null);
                return;
            }

            const matches: FindMatch[] = [];

            const colIndex = orderedParamColumns.findIndex((col) => col.field === findColumn);

            if (colIndex === -1) {
                console.log("Column index not found:", findColumn);
                return;
            }

            const column = orderedParamColumns[colIndex];

            if (!column) {
                console.log("Column not found:", findColumn);
                return;
            }

            for (let row = 0; row < filteredIndices.length; row++) {
                const dataIndex = filteredIndices[row];

                // orderedParamColumns.forEach((column, colIndex) => {
                const value = dataStore.getRowText(dataIndex, column.field);

                if (value) {
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
                // })
            }

            console.log(`Found ${matches.length} matches`);
            setFindMatches(matches);
            setMatchCount(matches.length);

            if (matches.length > 0) {
                setCurrentMatchIndex(0);
                const grid = gridRef.current?.slickGrid;
                if (grid) {
                    const match = matches[0];
                    // Set flag before navigation
                    isSelectingRef.current = true;

                    // grid.scrollRowIntoView(match.rowIndex, false);
                    // grid.scrollCellIntoView(match.rowIndex, match.columnIndex, false);
                    // grid.setActiveCell(match.rowIndex, match.columnIndex);

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
        [filteredIndices, orderedParamColumns, dataStore, config.include_index, findColumn],
    );

    const handleFindNext = useCallback(() => {
        // Sanity checks
        if (findMatches.length === 0) return;
        if (currentMatchIndex + 1 >= findMatches.length) return;

        const nextIndex = currentMatchIndex + 1;
        setCurrentMatchIndex(nextIndex);
        const grid = gridRef.current?.slickGrid;
        const match = findMatches[nextIndex];

        if (grid && match) {
            // Set flag before navigation
            isSelectingRef.current = true;

            grid.gotoCell(match.rowIndex, match.columnIndex, false);

            // Reset flag after delay
            setTimeout(() => {
                isSelectingRef.current = false;
            }, 100);
        }

        console.log(`Match ${nextIndex + 1} of ${findMatches.length}`);
    }, [findMatches, currentMatchIndex]);

    const handleFindPrev = useCallback(() => {
        // Sanity checks
        if (findMatches.length === 0) return;
        if (currentMatchIndex === 0) return;

        const prevIndex = (currentMatchIndex - 1) % findMatches.length;
        setCurrentMatchIndex(prevIndex);

        const grid = gridRef.current?.slickGrid;
        const match = findMatches[prevIndex];

        if (grid && match) {
            // Set flag before navigation
            isSelectingRef.current = true;

            grid.gotoCell(match.rowIndex, match.columnIndex, false);

            // Reset flag after delay
            setTimeout(() => {
                isSelectingRef.current = false;
            }, 100);
        }

        console.log(`Match ${prevIndex + 1} of ${findMatches.length}`);
    }, [findMatches, currentMatchIndex]);

    const options: GridOption = useMemo(
        () => ({
            gridWidth: "100%",
            gridHeight: "600px",
            // gridHeight: "100%",
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
            // editable: true,
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

    const onClose = useCallback(() => {
        setIsFindReplaceOpen(false);
        setFindColumn(null);
        setFindMatches([]);
        setMatchCount(null);
        setCurrentMatchIndex(-1);
    }, []);

    return (
        <div id={`react-grid-${chartId}`} className="absolute w-[100%] h-[100%] slickgrid-react-container">
            <SlickgridReact
                gridId={`table-${chartId}`}
                columns={columnDefs}
                dataset={[]}
                options={options}
                onReactGridCreated={handleGridCreated}
            />

            <FindAndReplaceDialog
                open={isFindReplaceOpen}
                onClose={onClose}
                handleFind={handleFind}
                handleFindPrev={handleFindPrev}
                handleFindNext={handleFindNext}
                handleReplace={(find: string, replace: string) => console.log("Replace: ", find, replace)}
                handleReplaceAll={(find: string, replace: string) => console.log("Replace All: ", find, replace)}
                foundMatches={matchCount}
                columnName={findColumn}
                // todo: might take a look at this to clean up
                disableFindPrev={findMatches.length === 0 || (findMatches.length > 0 && currentMatchIndex <= 0)}
                disableFindNext={
                    findMatches.length === 0 || (findMatches.length > 0 && currentMatchIndex + 1 >= findMatches.length)
                }
            />
        </div>
    );
});

export default TableChartReactComponent;
