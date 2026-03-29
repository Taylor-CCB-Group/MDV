import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { TableChartReact, TableChartReactConfig } from "../components/TableChartReactWrapper";
import { useChart, useDataStore } from "../context";
import { useChartID, useConfig, useOrderedParamColumns, useTheme } from "../hooks";
import { useHighlightedIndices } from "../selectionHooks";
import { type Column, Editors, type GridOption, type SlickgridReactInstance, SlickEventHandler } from "slickgrid-react";
import SlickGridDataProvider from "../utils/SlickGridDataProvider";
import { action, runInAction } from "mobx";
import useSortedFilteredIndices from "./useSortedFilteredIndices";

import { InputEditor } from "slickgrid-react";
import { DataModel } from "@/table/DataModel";
import type { FeedbackAlert } from "../components/FeedbackAlertComponent";
import type { AddColumnParams } from "../components/AddTableColumnDialog";
import type { BulkEditAction } from "../components/BulkEditColumnDialog";

/**
 * Text editor that sets the HTML input maxLength so the user cannot type
 * more than the allowed characters. Used by unique columns as of now.
 */
export class TextEditorWithMaxLength extends InputEditor {
    init(): void {
        super.init();
        const maxLen = this.columnEditor?.maxLength;
        const input = this.editorDomElement as HTMLInputElement | undefined;
        if (maxLen != null && input) {
            // Limit the characters to max length
            input.setAttribute("maxLength", String(maxLen));
        }
    }
}
/**
 * Hook for managing SlickGrid React component state.
 * 
 * State Management Contract:
 * 
 * 1. Config Object (MobX Observable):
 *    - The config object is deeply observable via MobX (made observable in BaseChart constructor)
 *    - All config properties should exist before makeAutoObservable is called (handled in adaptConfig)
 *    - Config properties are the source of truth for persisted state (column_widths, sort, order, etc.)
 * 
 * 2. Column Widths State Management:
 *    - config.column_widths is used ONLY for initial column widths when creating columnDefs
 *    - After grid is created, grid manages widths internally
 *    - columnDefs does NOT depend on config.column_widths (removed from dependencies)
 *    - getConfig() serializes current grid widths to config for persistence
 *    - No runtime synchronization needed - grid is authoritative during runtime
 * 
 * 3. Sort State Management:
 *    - suppressSortSyncRef prevents feedback loops
 *    - Config.sort is the source of truth, synced bidirectionally with grid visual state
 * 
 * 4. Event Handler Lifecycle:
 *    - Event handlers are attached in a useEffect that depends on [config, gridInstance]
 *    - Handlers update config properties for persistence
 *    - Suppression flags (only for sort) prevent reactive recalculation during operations
 */
const useSlickGridReact = () => {
    // Hooks
    const config = useConfig<TableChartReactConfig>();
    const dataStore = useDataStore();
    const chartId = useChartID();
    const chart = useChart<TableChartReactConfig, TableChartReact>();
    const orderedParamColumns = useOrderedParamColumns<TableChartReactConfig>();
    const sortedFilteredIndices = useSortedFilteredIndices();
    const highlightedIndices = useHighlightedIndices();
    const theme = useTheme();

    // States
    const [isFindReplaceOpen, setIsFindReplaceOpen] = useState(false);
    const [searchColumn, setSearchColumn] = useState<string | null>(null);
    const [feedbackAlert, setFeedbackAlert] = useState<FeedbackAlert>(null);
    const [isAddColumnDialogOpen, setIsAddColumnDialogOpen] = useState(false);
    const [isBulkEditDialogOpen, setIsBulkEditDialogOpen] = useState(false);
    const [bulkEditColumn, setBulkEditColumn] = useState<string | null>(null);

    // Refs
    const sortedFilteredIndicesRef = useRef(sortedFilteredIndices); // Holds latest value of indices when event handlers are called
    const orderedParamColumnsRef = useRef(orderedParamColumns); //  Holds latest value of columns when event handlers are called
    const selectionSourceRef = useRef<'user' | 'programmatic' | null>(null); // Track source of selection changes
    const gridRef = useRef<SlickgridReactInstance | null>(null);
    const suppressSortSyncRef = useRef(false); // Flag to prevent feedback loops during sort sync
    const cleanupRef = useRef<(() => void) | null>(null); // Cleanup event handlers

    useEffect(() => {
        sortedFilteredIndicesRef.current = sortedFilteredIndices;
        orderedParamColumnsRef.current = orderedParamColumns;
    }, [sortedFilteredIndices, orderedParamColumns]);

    useEffect(() => {
        chart.setAddColumnDialogOpener(() => setIsAddColumnDialogOpen(true));
        return () => {
            chart.setAddColumnDialogOpener(undefined);
        };
    }, [chart]);

    // Extract initial widths from config once (non-observable)
    // This avoids rules-of-hooks issues and prevents columnDefs from reacting to config.column_widths changes
    // We read config.column_widths only once on mount to get initial values
    const initialWidths = useMemo(() => {
        return config.column_widths ? { ...config.column_widths } : {};
        // we don't expect this to re-run, nothing should be touching `column_widths` at runtime.
        // keeping in dependency array to satisfy lint.
    }, [config.column_widths]);

    /**
     * Column definitions for SlickGrid.
     * 
     * State Management Contract:
     * - Uses initialWidths (extracted from config.column_widths on mount) ONLY for initial column widths
     * - After grid is created, grid manages widths internally
     * - getConfig() serializes current grid widths to config for persistence
     * - No runtime synchronization - grid is authoritative during runtime
     */
    const columnDefs = useMemo<Column[]>(() => {
        const cols: Column[] = [];

        if (config.include_index) {
            const indexWidth = initialWidths["__index__"] ?? 100;
            cols.push({
                id: "__index__",
                field: "__index__",
                name: "index",
                sortable: true,
                width: indexWidth,
                reorderable: false,
            });
        }

        for (const col of orderedParamColumns) {
            const isColumnEditable = col?.editable ?? false;
            const colWidth = initialWidths[col.field] ?? 100;
            cols.push({
                id: col.field,
                field: col.field,
                name: isColumnEditable
                ? `<span class="mdv-col-editable-icon mdi mdi-square-edit-outline" title="This column is editable, click on cell to edit" aria-hidden="true"></span><span class="mdv-col-name-text">${col.name}</span>`
                : col.name,
                sortable: true,
                width: colWidth,
                editor: isColumnEditable
                    ? col.datatype === "unique" && col.stringLength
                        // using a custom text editor to limit the characters to string length
                        //! we might need to think of increasing the string length in this case
                        ? { model: TextEditorWithMaxLength, maxLength: col.stringLength }
                        : { model: Editors.text }
                    : null,
                cssClass: isColumnEditable ? "mdv-editable-cell" : "",
                header: {
                    menu: {
                        commandItems: [
                            {
                                command: "find-replace",
                                title: isColumnEditable ? "Find & Replace" : "Find",
                                iconCssClass: "mdi mdi-magnify",
                            },
                            ...(isColumnEditable
                                ? [
                                    {
                                        command: "bulk-edit",
                                        title: "Bulk Edit",
                                        iconCssClass: "mdi mdi-table-edit",
                                    },
                                ]
                                : []),
                        ],
                    },
                },
            });
        }

        return cols;
    }, [config.include_index, orderedParamColumns, initialWidths]);

    /**
     * A new data provider gets created whenever any of the dependencies change which
     * will trigger an update to the grid
     * 
     * Creation of a new data provider isn't expensive according to LLM
     */
    const dataProvider = useMemo(() => {
        return new SlickGridDataProvider(orderedParamColumns, sortedFilteredIndices, config.include_index);
    }, [orderedParamColumns, sortedFilteredIndices, config.include_index]);

    /**
     * This useEffect is called whenever there is an update in dataProvider
     * (sorting, filtering, editing, etc)
     */
    useEffect(() => {
        const grid = gridRef.current?.slickGrid;
        if (grid && dataProvider) {
            // Skip update if selectionSourceRef is not empty
            if (selectionSourceRef.current !== null) {
                console.log("Skipping grid update during user selection");
                return;
            }
            grid.setData(dataProvider, true);
            grid.invalidate();
            console.log("Grid updated");
        }
    }, [dataProvider]);

    /**
     * This function handles grid creation, attaching data provider to the grid
     * and attaching all the event handlers
     * 
     * Only runs when onReactGridCreated event is called by grid
     */
    const handleGridCreated = useCallback(
        (e: CustomEvent<SlickgridReactInstance>) => {
            console.log("Grid created");
            gridRef.current = e.detail;
            
            // Store gridRef in chart instance for getConfig() access
            chart.setGridRef(gridRef);

            const grid = e.detail.slickGrid;
            if (grid && dataProvider) {
                grid.setData(dataProvider, true);
                grid.render();
                
                // Apply initial sort if config.sort is set
                // we need this because the useEffect for sort sync won't run if grid is null
                // and it won't rerun if we don't initially set the config.sort
                if (config.sort) {
                    suppressSortSyncRef.current = true;
                    grid.setSortColumn(config.sort.columnId, config.sort.ascending);
                    suppressSortSyncRef.current = false;
                }

                attachEventHandlers(e.detail);
            }
        },
        [dataProvider, chart, config.sort],
    );

    /**
     * Attach the event handlers to the grid
     * 
     * - Making use of a single SlickEventHandler instance for subscribing to the event to 
     * easily unsubscribe to all events except the pubService events
     * - Making use of sortedFilteredIndicesRef to get the current value when the handler is called
     * - selectionSourceRef keeps track of the selectedRows
     */
    const attachEventHandlers = useCallback((gridInstance: SlickgridReactInstance) => {
        const grid = gridInstance?.slickGrid;
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
                dataStore.dataHighlighted(indices, chart);
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
                const currentSort = {columnId, ascending: sortAsc};

                const currentSortStr = JSON.stringify(currentSort);
                const configSortStr = JSON.stringify(config.sort);
                
                const isSame = currentSortStr === configSortStr;
                console.log("Sort event:", columnId, sortAsc ? "asc" : "desc");
                // As far as the types go... I think if the `sortAsc` is undefined, then that means `config.sort` should be undefined.
                // if not, then the type of config.sort.ascending should be optional.
                // casting `as bool` just means that we make the types lie.
                // It may seem like we don't reach this case because of `"sortAsc" in args` above - but technically,
                // this gets into fiddly "key-optional" vs "value-optional" distinction.
                // args.sortAsc could have a value of undefined - that's different from the key not being in the object in a meaningful way.
                // it shouldn't be treated the same as a false value, or typed as though it was and passed further down...
                if (!isSame) { // Only update if there is a change
                    if (sortAsc === undefined) {
                        config.sort = undefined;
                    } else {
                        config.sort = { columnId, ascending: sortAsc };
                    }
                }
            }
        }));

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
                } else if (command === "bulk-edit") {
                    setBulkEditColumn(column.field);
                    setIsBulkEditDialogOpen(true);
                }
            },
        ));

        const gridMenuSubscription = pubSub.subscribe("onGridMenuCommand", action(({ command }: { command: string }) => {
            if (command === "clear-sorting") {
                config.sort = undefined;
            }
        }));

        cleanupRef.current = () => {
            slickEventHandler.unsubscribeAll();
            // the typing isn't great for these - but we believe as of writing that these will exist, 
            // and unsubscribe() will be called on both.
            if (!headerMenuSubscription) {
                console.warn("no headerMenuSubscription to unsubscribe from... minor leak here.");
            }
            headerMenuSubscription?.unsubscribe?.();
            gridMenuSubscription?.unsubscribe?.();
        };
    }, [config, chart, dataStore]);

    useEffect(() => {
        // Cleanup the event handlers on unmount
        return () => {
            if (cleanupRef.current) {
                cleanupRef.current();
            }
        }
    }, []);

    /**
     * Sync config.sort → grid visual state
     * 
     * If the sort is changed externally or internally through grid, this synchronizes 
     * the config.sort
     * suppressSortSyncRef avoids multiple rerenders when the sortColumns change
     */
    useEffect(() => {
        const grid = gridRef.current?.slickGrid;
        if (!grid) return;
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
    }, [config.sort]);

    /**
     * Handle external highlighted data
     * 
     * When highlighting changes externally, update the grid's visual selection
     * Internal grid's highlighting is handled by the event handler
     * 
     */
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

    const isColumnEditable = useMemo(() => {
        const column = orderedParamColumns.find((col) => col.field === searchColumn);
        return column?.editable ?? false;
    }, [searchColumn, orderedParamColumns]);

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

    const closeAddColumnDialog = useCallback(() => {
        setIsAddColumnDialogOpen(false);
    }, []);

    const closeBulkEditDialog = useCallback(() => {
        setIsBulkEditDialogOpen(false);
        setBulkEditColumn(null);
    }, []);

    // Columns to be displayed for cloning
    const cloneableColumns = useMemo(() => {
        return orderedParamColumns
            .map((column) => ({
                field: column.field,
                name: column.name,
                datatype: column.datatype,
                stringLength: column.stringLength,
                delimiter: column.delimiter,
            }))
            .sort((a, b) => a.name.localeCompare(b.name));
    }, [orderedParamColumns]);

    const addColumnDefaultPosition = useMemo(() => {
        return config.include_index ? 2 : Math.min(2, orderedParamColumns.length + 1);
    }, [config.include_index, orderedParamColumns.length]);

    const handleAddColumn = useCallback(({
        name,
        datatype,
        cloneColumn,
        position,
        stringLength,
        delimiter,
    }: AddColumnParams) => {
        const trimmedName = name.trim();

        if (!trimmedName) {
            setFeedbackAlert({
                type: "warning",
                title: "Add Column Warning",
                message: "Column name is required.",
            });
            return;
        }

        if (dataStore.columnIndex[trimmedName]) {
            setFeedbackAlert({
                type: "error",
                title: "Add Column Error",
                message: `Column ${trimmedName} already exists`,
            });
            return;
        }

        try {
            const dataModel = new DataModel(dataStore, { autoupdate: false });
            // Create a new column
            dataModel.createColumn({
                name: trimmedName,
                datatype: datatype,
                cloneColumn,
                stringLength,
                delimiter,
            });

            // Get the position of the column to be inserted
            const orderedFields = orderedParamColumnsRef.current.map((column) => column.field);
            const displayedFields = config.include_index
                ? ["__index__", ...orderedFields]
                : [...orderedFields];
            const parsedPosition = position ?? displayedFields.length + 1;
            const insertionIndex = Math.min(
                Math.max(parsedPosition - 1, 0),
                displayedFields.length,
            );

            // Insert the column in the position and update the order and param fields
            displayedFields.splice(insertionIndex, 0, trimmedName);
            const nextFields = displayedFields.filter((field) => field !== "__index__");
            const nextOrder = Object.fromEntries(
                nextFields.map((field, index) => [field, index]),
            );

            runInAction(() => {
                config.param = nextFields;
                config.order = nextOrder;
            });

            // Notify datastore and rerender the grid
            dataStore.dataChanged([trimmedName]);
            gridRef.current?.slickGrid?.invalidate();
            setIsAddColumnDialogOpen(false);
        } catch (err) {
            const error =
                err instanceof Error ? err : new Error("Failed to add editable column");
            setFeedbackAlert({
                type: "error",
                title: "Add Column Error",
                message: error.message,
                stack: error.stack,
                metadata: {
                    columnName: trimmedName,
                    cloneColumn,
                    position,
                },
            });
        }
    }, [config, dataStore]);

    const handleBulkEdit = useCallback(({
        action,
        columnName,
        value,
    }: {
        action: BulkEditAction;
        columnName: string;
        value: string;
    }) => {
        try {
            const dataModel = new DataModel(dataStore, { autoupdate: false });
            const rowIndices = Array.from(sortedFilteredIndicesRef.current);

            if (action === "fill-all") {
                dataModel.fillColumn(columnName, value, rowIndices, false);
            } else {
                dataModel.fillColumn(columnName, value, rowIndices, true);
            }

            gridRef.current?.slickGrid?.invalidate();
            closeBulkEditDialog();
        } catch (err) {
            const error =
                err instanceof Error ? err : new Error("Failed to bulk edit column");
            setFeedbackAlert({
                type: "error",
                title: "Bulk Edit Error",
                message: error.message,
                stack: error.stack,
                metadata: {
                    action,
                    columnName,
                },
            });
        }
    }, [closeBulkEditDialog, config, dataStore]);

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
        feedbackAlert,
        setFeedbackAlert,
        isAddColumnDialogOpen,
        cloneableColumns,
        addColumnDefaultPosition,
        closeAddColumnDialog,
        handleAddColumn,
        isBulkEditDialogOpen,
        bulkEditColumn,
        closeBulkEditDialog,
        handleBulkEdit,
    };
};

export default useSlickGridReact;
