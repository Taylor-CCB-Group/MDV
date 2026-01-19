import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { TableChartReactConfig } from "../components/TableChartReactWrapper";
import { useChart, useDataStore } from "../context";
import { useChartID, useConfig, useOrderedParamColumns, useTheme } from "../hooks";
import { useHighlightedIndices } from "../selectionHooks";
import {
    type Column,
    Editors,
    type GridOption,
    type SlickgridReactInstance,
    SlickEventHandler,
} from "slickgrid-react";
import SlickGridDataProvider from "../utils/SlickGridDataProvider";
import { runInAction } from "mobx";
import useSortedFilteredIndices from "./useSortedFilteredIndices";

const useSlickGridReact = () => {
    // Hooks
    const config = useConfig<TableChartReactConfig>();
    const dataStore = useDataStore();
    const chartId = useChartID();
    const chart = useChart<TableChartReactConfig>();
    const orderedParamColumns = useOrderedParamColumns<TableChartReactConfig>();
    const sortedFilteredIndices = useSortedFilteredIndices();
    const highlightedIndices = useHighlightedIndices();
    const theme = useTheme();

    // States
    const [isFindReplaceOpen, setIsFindReplaceOpen] = useState(false);
    const [searchColumn, setSearchColumn] = useState<string | null>(null);
    const [gridInstance, setGridInstance] = useState<SlickgridReactInstance | null>(null);

    // Refs
    const sortedFilteredIndicesRef = useRef(sortedFilteredIndices);
    const orderedParamColumnsRef = useRef(orderedParamColumns);
    const dataStoreRef = useRef(dataStore);
    const chartRef = useRef(chart);
    const isSelectingRef = useRef(false); // Flag to prevent grid update during selection
    const gridRef = useRef<SlickgridReactInstance | null>(null);

    useEffect(() => {
        sortedFilteredIndicesRef.current = sortedFilteredIndices;
        orderedParamColumnsRef.current = orderedParamColumns;
        dataStoreRef.current = dataStore;
        chartRef.current = chart;
    }, [sortedFilteredIndices, dataStore, chart, orderedParamColumns]);

    const columnDefs = useMemo<Column[]>(() => {
        const cols: Column[] = [];

        if (config.include_index) {
            cols.push({
                id: "__index__",
                field: "__index__",
                name: "index",
                sortable: true,
                width: config.column_widths?.["__index__"] ?? 100,
            });
        }

        for (const col of orderedParamColumns) {
            const isColumnEditable = col?.editable ?? false;
            cols.push({
                id: col.field,
                field: col.field,
                name: col.name,
                sortable: true,
                width: config.column_widths?.[col.field] ?? 100,
                editor: isColumnEditable ? { model: Editors.text } : null,
                cssClass: isColumnEditable ? "mdv-editable-cell" : "",
                header: {
                    menu: {
                        commandItems: [
                            {
                                command: "find-replace",
                                title: isColumnEditable ? "Find & Replace" : "Find",
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
        return new SlickGridDataProvider(orderedParamColumns, sortedFilteredIndices, config.include_index);
    }, [orderedParamColumns, sortedFilteredIndices, config.include_index]);

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

            setGridInstance(e.detail);
            const grid = e.detail.slickGrid;
            if (grid && dataProvider) {
                grid.setData(dataProvider, true);
                grid.render();
            }
        },
        [dataProvider],
    );

    useEffect(() => {
        // Using grid state to attach the handlers after the grid is created
        if (!gridInstance) return;
        const grid = gridInstance.slickGrid;
        if (!grid) return;
        const slickEventHandler = new SlickEventHandler();

        slickEventHandler.subscribe(grid.onSelectedRowsChanged, (_e, args) => {
            const selectedRows = args.rows;

            if (selectedRows.length > 0) {
                const indices = selectedRows.map((row) => sortedFilteredIndicesRef.current[row]);

                // Set flag to prevent grid update during selection
                isSelectingRef.current = true;

                dataStoreRef.current.dataHighlighted(indices, chartRef.current);

                // Reset flag after reactions settle
                setTimeout(() => {
                    isSelectingRef.current = false;
                }, 100);
            }
        });

        slickEventHandler.subscribe(grid.onSort, (_e, args) => {
            if ("sortCol" in args && args.sortCol && "sortAsc" in args) {
                const columnId = args.sortCol.field;
                const sortAsc = args.sortAsc;
                console.log("Sort event:", columnId, sortAsc ? "asc" : "desc");
                runInAction(() => {
                    // As far as the types go... I think if the `sortAsc` is undefined, then that means `config.sort` should be undefined.
                    // if not, then the type of config.sort.ascending should be optional.
                    // casting `as bool` just means that we make the types lie.
                    // It may seem like we don't reach this case because of `"sortAsc" in args` above - but technically,
                    // this gets into fiddly "key-optional" vs "value-optional" distinction.
                    // args.sortAsc could have a value of undefined - that's different from the key not being in the object in a meaningful way.
                    // it shouldn't be treated the same as a false value, or typed as though it was and passed further down...
                    if (sortAsc === undefined) {
                        config.sort = undefined;
                    } else {
                        config.sort = { columnId, ascending: sortAsc };
                    }
                });
            }
        });

        const pubSub = grid.getPubSubService();
        if (!pubSub) {
            // strong assertion acts as a type-guard
            throw new Error("SlickGrid PubSubService undefined - this should be unreachable.");
        }
        // don't think we can use our slickEventHandler with pubSub
        const headerMenuSubscription = pubSub.subscribe("onHeaderMenuCommand", (event: { column: Column; command: string; }) => {
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

        const gridMenuSubscription = pubSub.subscribe("onGridMenuCommand", ({ command }: { command: string; }) => {
            if (command === "clear-sorting") {
                console.log("Clear All sort");
                runInAction(() => {
                    config.sort = undefined;
                });
            }
        });
        
        
        slickEventHandler.subscribe(grid.onColumnsResized, (_e, args) => {
            console.log("resize handler");
            const columns = args.grid.getColumns();
            runInAction(() => {
                // Create a new reference of config.columnWidths
                const columnWidths: Record<string, number> = config.column_widths ? { ...config.column_widths } : {};

                columns.forEach((col: Column) => {
                    if (col.width && col.width !== 100) {
                        columnWidths[col.field] = col.width;
                    } else if (columnWidths[col.field]) {
                        // Remove if reset to default
                        delete columnWidths[col.field];
                    }
                });

                // Change the reference of config.order for react to detect and update
                config.column_widths = columnWidths;
            });
        });

        slickEventHandler.subscribe(grid.onColumnsReordered, (_e, args) => {
            console.log("reorder handler");
            const impactedColumns = args.impactedColumns;
            runInAction(() => {
                // Create a new reference of config.order
                const newOrder = config.order ? { ...config.order } : {};
                const cols = impactedColumns.filter((col) => col.field !== "__index__");
                cols.forEach((col, index) => {
                    newOrder[col.field] = index;
                });

                // Change the reference of config.order for react to detect and update
                config.order = newOrder;
            });
        });

        return () => {
            slickEventHandler.unsubscribeAll();
            headerMenuSubscription.unsubscribe?.();
            gridMenuSubscription.unsubscribe?.();
        };
    }, [config, gridInstance]);

    // Handle highlighted data
    useEffect(() => {
        const grid = gridRef.current?.slickGrid;
        if (!grid) return;

        if (!dataProvider || !sortedFilteredIndices || sortedFilteredIndices.length === 0) {
            return;
        }

        try {

        if (highlightedIndices.length === 0) {
            // Only reset if the data provider is initialized and sorted indices have values
            // otherwise we will be messing with the initialization of the grid
            isSelectingRef.current = true;
            grid.setSelectedRows([]);
            setTimeout(() => {
                isSelectingRef.current = false;
            }, 100);
            return;
        }

        const filteredSet = new Set(sortedFilteredIndicesRef.current);

        // Filter the highlightedIndices by checking if the sortedFilteredIndices have those indices
        const validIndices = highlightedIndices.filter(i => filteredSet.has(i));

        if (validIndices.length === 0) {
            isSelectingRef.current = true;
            grid.setSelectedRows([]);
            setTimeout(() => {
                isSelectingRef.current = false;
            }, 100);
            return;
        }

        const positions: number[] = [];

        for (const index of highlightedIndices) {
            // Get the position of the index from sorted indices
            const pos = sortedFilteredIndicesRef.current.indexOf(index);
            if (pos !== -1) positions.push(pos);
        }

        if (positions.length === 0) return;

        isSelectingRef.current = true;

        // Navigate to the first row
        grid.scrollRowIntoView(positions[0], false);
        // Set the selected rows in the grid
        grid.setSelectedRows(positions);

        setTimeout(() => {
            isSelectingRef.current = false;
        }, 100)
    } catch (err) {
        console.error("Error highlighting the rows in the table", err);
    }
    }, [highlightedIndices, dataProvider, sortedFilteredIndices]);

    const options: GridOption = useMemo(
        () => ({
            gridWidth: "100%",
            gridHeight: "600px",
            darkMode: theme === "dark",
            autoFitColumnsOnFirstLoad: false, // To avoid the columns to take less width
            enableAutoSizeColumns: false, 
            rowHeight: 25,  
            defaultColumnWidth: 100,
            enableColumnPicker: false, // Disable right click on column header
            enableSorting: true,
            multiColumnSort: false,
            enableGridMenu: false, // Disabled as it's interfering with the last column
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
            resizeByContentOptions: {
                // hack - seems somewhat closer to getting the right size here?
                defaultRatioForStringType: 1.1,
            },
            editable: true,
            enableCellNavigation: true,
            enableExcelCopyBuffer: true,
            multiSelect: true,
            enableRowSelection: true,
            rowSelectionOptions: {
                selectActiveRow: true,
            },
            headerMenu: {
                hideColumnHideCommand: true, // Disabled to avoid messing with the column order
            },
        } satisfies GridOption), // helps with autocomplete for options
        [chartId, theme],
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
        sortedFilteredIndices,
        isFindReplaceOpen,
        searchColumn,
        setIsFindReplaceOpen,
        setSearchColumn,
        sortedFilteredIndicesRef,
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
