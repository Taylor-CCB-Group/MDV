import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { TableChartReactConfig } from "../components/TableChartReactWrapper";
import { useChart, useDataStore } from "../context";
import { useChartID, useConfig, useOrderedParamColumns, useTheme } from "../hooks";
import { useHighlightedIndices } from "../selectionHooks";
import { type Column, Editors, type GridOption, type SlickgridReactInstance, SlickEventHandler } from "slickgrid-react";
import SlickGridDataProvider from "../utils/SlickGridDataProvider";
import { action } from "mobx";
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
    const selectionSourceRef = useRef<'user' | 'programmatic' | null>(null); // Track source of selection changes
    const gridRef = useRef<SlickgridReactInstance | null>(null);
    const suppressSortSyncRef = useRef(false); // Flag to prevent feedback loops during sort sync

    useEffect(() => {
        sortedFilteredIndicesRef.current = sortedFilteredIndices;
        orderedParamColumnsRef.current = orderedParamColumns;
        dataStoreRef.current = dataStore;
        chartRef.current = chart;
    }, [sortedFilteredIndices, dataStore, chart, orderedParamColumns]);

    const columnDefs = useMemo<Column[]>(() => {
        const cols: Column[] = [];

        // Get current column widths from grid if it exists, to preserve autosized widths
        // that might not yet be saved to config.column_widths
        let currentGridWidths: Record<string, number> | null = null;
        const grid = gridRef.current?.slickGrid;
        if (grid) {
            try {
                const gridColumns = grid.getColumns();
                currentGridWidths = {};
                for (const col of gridColumns) {
                    if (col.width && col.width !== 100) {
                        currentGridWidths[col.field] = col.width;
                    }
                }
            } catch (e) {
                // Grid might not be fully initialized, ignore
                console.debug("Could not read current grid widths:", e);
            }
        }

        if (config.include_index) {
            const indexWidth = currentGridWidths?.["__index__"] ?? config.column_widths?.["__index__"] ?? 100;
            cols.push({
                id: "__index__",
                field: "__index__",
                name: "index",
                sortable: true,
                width: indexWidth,
            });
        }

        for (const col of orderedParamColumns) {
            const isColumnEditable = col?.editable ?? false;
            // Prefer current grid width (if exists) over saved config width, to preserve autosized widths
            const colWidth = currentGridWidths?.[col.field] ?? config.column_widths?.[col.field] ?? 100;
            cols.push({
                id: col.field,
                field: col.field,
                name: col.name,
                sortable: true,
                width: colWidth,
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
            // Skip update if user is actively selecting
            if (selectionSourceRef.current === 'user') {
                console.log("Skipping grid update during user selection");
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
                
                // Apply initial sort if config.sort is set
                if (config.sort) {
                    suppressSortSyncRef.current = true;
                    grid.setSortColumn(config.sort.columnId, config.sort.ascending);
                    suppressSortSyncRef.current = false;
                }
            }
        },
        [dataProvider, config.sort],
    );

    useEffect(() => {
        // Using grid state to attach the handlers after the grid is created
        if (!gridInstance) return;
        const grid = gridInstance.slickGrid;
        if (!grid) return;
        const slickEventHandler = new SlickEventHandler();

        slickEventHandler.subscribe(grid.onSelectedRowsChanged, (_e, args) => {
            const selectedRows = args.rows;

            // Skip if we're programmatically updating
            if (selectionSourceRef.current === 'programmatic') {
                selectionSourceRef.current = null;
                return;
            }

            if (selectedRows.length > 0) {
                selectionSourceRef.current = 'user';
                const indices = selectedRows.map((row) => sortedFilteredIndicesRef.current[row]);
                dataStoreRef.current.dataHighlighted(indices, chartRef.current);
                // Reset immediately - the effect will handle any needed updates
                selectionSourceRef.current = null;
            }
        });

        slickEventHandler.subscribe(grid.onSort, action((_e, args) => {
            // Skip during programmatic sync to prevent feedback loops
            if (suppressSortSyncRef.current) return;
            
            if ("sortCol" in args && args.sortCol && "sortAsc" in args) {
                const columnId = args.sortCol.field;
                const sortAsc = args.sortAsc;
                console.log("Sort event:", columnId, sortAsc ? "asc" : "desc");
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
            }
        }));

        const pubSub = grid.getPubSubService();
        if (!pubSub) {
            // strong assertion acts as a type-guard
            throw new Error("SlickGrid PubSubService undefined - this should be unreachable.");
        }
        // don't think we can use our slickEventHandler with pubSub
        const headerMenuSubscription = pubSub.subscribe(
            "onHeaderMenuCommand",
            action((event: { column: Column; command: string }) => {
                const { column, command } = event;
                if (command === "clear-sort") {
                    config.sort = undefined;
                } else if (command === "sort-asc") {
                    config.sort = { columnId: column.field, ascending: true };
                } else if (command === "sort-desc") {
                    config.sort = { columnId: column.field, ascending: false };
                } else if (command === "find-replace") {
                    setSearchColumn(column.field);
                    setIsFindReplaceOpen(true);
                }
            },
        ));

        const gridMenuSubscription = pubSub.subscribe("onGridMenuCommand", action(({ command }: { command: string }) => {
            if (command === "clear-sorting") {
                config.sort = undefined;
            }
        }));

        // Subscribe to autosize event to capture autosized column widths.
        // This ensures autosized widths are saved immediately when autosizing occurs.
        slickEventHandler.subscribe(
            grid.onAutosizeColumns,
            action((_e, _args) => {
                const columns = grid.getColumns();
                const columnWidths: Record<string, number> = config.column_widths
                    ? { ...config.column_widths }
                    : {};

                for (const col of columns) {
                    if (col.width && col.width !== 100) {
                        // Save any non-default width (including autosized widths)
                        columnWidths[col.field] = col.width;
                    } else if (columnWidths[col.field]) {
                        // Remove if reset to default
                        delete columnWidths[col.field];
                    }
                }

                config.column_widths = columnWidths;
            }),
        );

        slickEventHandler.subscribe(
            grid.onColumnsResized,
            action((_e, args) => {
                // Save column widths when columns are resized (manually or via autosize)
                // This handler captures both manual resizes and autosizes that trigger this event
                const columns = args.grid.getColumns();
                const columnWidths: Record<string, number> = config.column_widths
                    ? { ...config.column_widths }
                    : {};

                for (const col of columns) {
                    if (col.width && col.width !== 100) {
                        columnWidths[col.field] = col.width;
                    } else if (columnWidths[col.field]) {
                        delete columnWidths[col.field];
                    }
                }

                config.column_widths = columnWidths;
            }),
        );

        slickEventHandler.subscribe(grid.onColumnsReordered, action((_e, args) => {
            const impactedColumns = args.impactedColumns;
            // Create a new reference of config.order
            const newOrder = config.order ? { ...config.order } : {};
            const cols = impactedColumns.filter((col) => col.field !== "__index__");
            cols.forEach((col, index) => {
                newOrder[col.field] = index;
            });

            // Change the reference of config.order for react to detect and update
            config.order = newOrder;
        }));

        return () => {
            slickEventHandler.unsubscribeAll();
            // the typing isn't great for these - but we believe as of writing that these will exist, 
            // and unsubscribe() will be called on both.
            if (!headerMenuSubscription) {
                console.warn("no headerMenuSubscription to unsubscribe from... minor leak here.");
            }
            headerMenuSubscription?.unsubscribe?.();
            gridMenuSubscription?.unsubscribe?.();
        };
    }, [config, gridInstance]);

    // Sync config.sort → grid visual state
    useEffect(() => {
        const grid = gridRef.current?.slickGrid;
        if (!grid || !gridInstance) return;
        if (suppressSortSyncRef.current) return;
        
        const currentSortCols = grid.getSortColumns();
        const currentSort = currentSortCols.length > 0 
            ? { columnId: currentSortCols[0].columnId, ascending: currentSortCols[0].sortAsc }
            : undefined;
        
        // Only update if different
        const configStr = JSON.stringify(config.sort);
        const currentStr = JSON.stringify(currentSort);
        
        if (configStr !== currentStr) {
            suppressSortSyncRef.current = true;
            try {
                if (config.sort) {
                    grid.setSortColumn(config.sort.columnId, config.sort.ascending);
                } else {
                    grid.setSortColumns([]);
                }
            } finally {
                suppressSortSyncRef.current = false;
            }
        }
    }, [config.sort, gridInstance]);

    // Handle highlighted data
    useEffect(() => {
        const grid = gridRef.current?.slickGrid;
        if (!grid) return;

        if (!dataProvider || !sortedFilteredIndices || sortedFilteredIndices.length === 0) {
            return;
        }

        try {
            // Get current selection from grid
            const currentSelection = grid.getSelectedRows();
            const currentSelectionSorted = [...currentSelection].sort();
            
            // Calculate desired selection
            const filteredSet = new Set(sortedFilteredIndices);
            const validIndices = highlightedIndices.filter((i) => filteredSet.has(i));
            
            if (validIndices.length === 0) {
                // Only clear if there's actually a selection and it wasn't user-initiated
                if (currentSelection.length > 0 && selectionSourceRef.current !== 'user') {
                    selectionSourceRef.current = 'programmatic';
                    grid.setSelectedRows([]);
                    requestAnimationFrame(() => {
                        selectionSourceRef.current = null;
                    });
                }
                return;
            }

            const positions = validIndices
                .map(index => sortedFilteredIndices.indexOf(index))
                .filter(pos => pos !== -1);

            if (positions.length === 0) return;
            
            const desiredSelectionSorted = [...positions].sort();
            
            // Only update if different
            if (JSON.stringify(currentSelectionSorted) !== JSON.stringify(desiredSelectionSorted)) {
                selectionSourceRef.current = 'programmatic';
                grid.scrollRowIntoView(positions[0], false);
                grid.setSelectedRows(positions);
                requestAnimationFrame(() => {
                    selectionSourceRef.current = null;
                });
            }
        } catch (err) {
            console.error("Error highlighting the rows in the table", err);
        }
    }, [highlightedIndices, dataProvider, sortedFilteredIndices]);

    const options: GridOption = useMemo(
        () =>
            ({
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
            }) satisfies GridOption, // helps with autocomplete for options
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
        selectionSourceRef,
        gridRef,
        options,
        columnDefs,
        handleGridCreated,
        isColumnEditable,
        onDialogClose,
    };
};

export default useSlickGridReact;
