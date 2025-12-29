import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { TableChartReactConfig } from "../components/TableChartReactWrapper";
import { useChart, useDataStore } from "../context";
import { useChartID, useConfig, useOrderedParamColumns } from "../hooks";
import { useHighlightedIndex } from "../selectionHooks";
import useSortedIndices from "./useSortedIndices";
import {
    type Column,
    Editors,
    type GridOption,
    type OnColumnsReorderedEventArgs,
    type OnColumnsResizedEventArgs,
    type SlickgridReactInstance,
} from "slickgrid-react";
import SlickGridDataProvider from "../utils/SlickGridDataProvider";
import { runInAction } from "mobx";

const useSlickGridReact = () => {
    // Hooks
    const config = useConfig<TableChartReactConfig>();
    const dataStore = useDataStore();
    const chartId = useChartID();
    const chart = useChart<TableChartReactConfig>();
    const orderedParamColumns = useOrderedParamColumns<TableChartReactConfig>();
    const sortedIndices = useSortedIndices();
    const highlightedIndex = useHighlightedIndex();

    // States
    const [isFindReplaceOpen, setIsFindReplaceOpen] = useState(false);
    const [searchColumn, setSearchColumn] = useState<string | null>(null);

    // Refs
    const sortedIndicesRef = useRef(sortedIndices);
    const orderedParamColumnsRef = useRef(orderedParamColumns);
    const dataStoreRef = useRef(dataStore);
    const chartRef = useRef(chart);
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
                name: "index",
                sortable: true,
                minWidth: config.column_widths?.["__index__"] || 100,
            });
        }

        for (const col of orderedParamColumns) {
            const isColumnEditable = col?.editable ?? false;
            cols.push({
                id: col.field,
                field: col.field,
                name: col.name,
                sortable: true,
                minWidth: config.column_widths?.[col.field] || 100,
                editor: isColumnEditable ? { model: Editors.text } : null,
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
        return new SlickGridDataProvider(orderedParamColumns, sortedIndices, config.include_index);
    }, [orderedParamColumns, sortedIndices, config.include_index]);

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

        const resizeHandler: any = grid.onColumnsResized?.subscribe((_e, args: OnColumnsResizedEventArgs) => {
            console.log("resize handler");
            const columns = args.grid.getColumns();
            runInAction(() => {
                if (!config.column_widths) {
                    config.column_widths = {};
                }
                // reassigning config.column_widths to local variable to avoid undefined type error
                const columnWidths = config.column_widths;
                columns.forEach((col: Column) => {
                    if (col.width && col.width !== 100) {
                        columnWidths[col.field] = col.width;
                    } else if (columnWidths[col.field]) {
                        // Remove if reset to default
                        delete columnWidths[col.field];
                    }
                });
            });
        });

        const reorderHandler: any = grid.onColumnsReordered?.subscribe((_e, args: OnColumnsReorderedEventArgs) => {
            console.log("reorder handler");
            const impactedColumns = args.impactedColumns;
            runInAction(() => {
                // Initialize order if it doesn't exist
                if (!config.order) {
                    config.order = {};
                }
                // reassigning config.order to local variable to avoid undefined type error
                const order = config.order;
                const cols = impactedColumns.filter((col) => col.field !== "__index__");
                cols.forEach((col, index) => {
                    order[col.field] = index;
                });
            });
        });

        return () => {
            selectionHandler?.unsubscribe();
            sortHandler?.unsubscribe();
            headerMenuHandler?.unsubscribe();
            gridMenuHandler?.unsubscribe();
            resizeHandler?.unsubscribe();
            reorderHandler?.unsubscribe();
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

    const onDialogClose = useCallback(() => {
        setIsFindReplaceOpen(false);
        setSearchColumn(null);
    }, []);

    return {
        config,
        dataStore,
        chartId,
        orderedParamColumns,
        sortedIndices,
        isFindReplaceOpen,
        searchColumn,
        setIsFindReplaceOpen,
        setSearchColumn,
        sortedIndicesRef,
        orderedParamColumnsRef,
        isSelectingRef,
        gridRef,
        options,
        columnDefs,
        handleGridCreated,
        isColumnEditable,
        onDialogClose,
    };
};

export default useSlickGridReact;
