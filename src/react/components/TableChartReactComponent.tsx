import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useMemo, useRef } from "react";
import { type Column, type GridOption, SlickgridReact, type SlickgridReactInstance } from "slickgrid-react";
import { useChartID, useConfig, useOrderedParamColumns } from "../hooks";
import type { TableChartReactConfig } from "./TableChartReactWrapper";
import SlickGridDataProvider from "../utils/SlickGridDataProvider";
import { useDataStore } from "../context";
import { runInAction } from "mobx";
import useSortedIndices from "../hooks/useSortedIndices";

// todo: Add find and replace and editing functionality
// todo: Fix the styling of the table to have the scroll bars visible
const TableChartReactComponent = observer(() => {
    const config = useConfig<TableChartReactConfig>();
    const dataStore = useDataStore();
    const chartId = useChartID();
    const orderedParamColumns = useOrderedParamColumns();
    const filteredIndices = useSortedIndices();

    const gridRef = useRef<SlickgridReactInstance | null>(null);

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
            });
        }

        return cols;
    }, [config.include_index, config.column_widths, orderedParamColumns]);

    const dataProvider = useMemo(() => {
        return new SlickGridDataProvider(dataStore, orderedParamColumns, filteredIndices, config.include_index);
    }, [dataStore, orderedParamColumns, filteredIndices, config.include_index]);

    const options: GridOption = useMemo(
        () => ({
            autoHeight: false,
            gridWidth: "100%",
            darkMode: window?.mdv?.chartManager?.theme === 'dark',
            enableSorting: true,
            multiColumnSort: false,
            alwaysShowVerticalScroll: true,
            alwaysAllowHorizontalScroll: true,
        }),
        [],
    );

    useEffect(() => {
        const grid = gridRef.current?.slickGrid;
        if (grid && dataProvider) {
            grid.setData(dataProvider, true);
            grid.render();
            console.log("Grid updated");
        }
    }, [dataProvider]);

    const handleGridCreated = useCallback((e: CustomEvent<SlickgridReactInstance>) => {
        console.log("Grid created");
        gridRef.current = e.detail;

        const grid = e.detail.slickGrid;
        if (grid && dataProvider) {
            grid.setData(dataProvider, true);
            grid.render();
        }
    }, [dataProvider]);

    useEffect(() => {
        const grid = gridRef.current?.slickGrid;
        if (!grid) return;

        const sortHandler: any = grid.onSort.subscribe((_e, args) => {
            if ("sortCol" in args && args.sortCol && "sortAsc" in args) {
                const columnId = args.sortCol.field as string;
                const sortAsc = args.sortAsc as boolean;
                console.log("Sort event:", columnId, sortAsc ? "asc" : "desc");
                runInAction(() => {
                    config.sort = { columnId, ascending: sortAsc };
                });
            }
        })

        const clearSortHandler = grid.getPubSubService()?.subscribe(
            "onHeaderMenuCommand",
            (event: { column: Column, command: string }) => {
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
                }
            }
        );

        const clearAllSortHandler = grid.getPubSubService()?.subscribe(
            "onGridMenuCommand",
            ({ command }: { command: string }) => {
                if (command === "clear-sorting") {
                    console.log("Clear All sort");
                    runInAction(() => {
                        config.sort = undefined;
                    });
                }
            }
        );

        return () => {
            sortHandler?.unsubscribe();
            clearSortHandler?.unsubscribe();
            clearAllSortHandler?.unsubscribe();
        }

    }, [config]);

    return (
        <div className="w-full h-full">
            <SlickgridReact
                gridId={`table-${chartId}`}
                columns={columnDefs}
                dataset={[]}
                options={options}
                onReactGridCreated={handleGridCreated}
            />
        </div>
    );
});

export default TableChartReactComponent;


