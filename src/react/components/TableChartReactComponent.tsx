import { observer } from "mobx-react-lite";
import { useCallback, useMemo, useState } from "react";
import { type Column, type GridOption, SlickgridReact } from "slickgrid-react";

// todo: Use the useOrderedParamColumns to get the columns
// Create a data provider probably to provide the data to the grid using index
// Add filtering, sorting, find and replace and editing functionality
// The logic would be quite similar to the SlickGridDemo code
// Modularize the code into different components

//* Trying to fix the styling of the grid by using mock data
const TableChartReactComponent = observer(() => {
    const [columns, setColumns] = useState<Column[]>([]);
    const [dataset, setDataset] = useState<any[]>([]);

    const gridOptions: GridOption = useMemo(
        () => ({
            // gridHeight: 450,
            // gridWidth: 800,
            // todo: fix dark mode to sync with the global mode
            darkMode: true,
            alwaysShowVerticalScroll: true,
            alwaysAllowHorizontalScroll: true,
        }),
        [],
    );

    const onGridCreated = useCallback(() => {
        const cols: Column[] = [
            { id: "id", name: "ID", field: "id", sortable: true, width: 70 },
            { id: "title", name: "Title", field: "title", sortable: true, width: 100 },
            { id: "duration", name: "Duration", field: "duration", sortable: true, width: 100 },
        ];
        const rows: any[] = [];
        for (let i = 0; i < 50; i++) {
            rows.push({
                id: i+1,
                title: `Task ${i+1}`,
                duration: `${i+1} days`,
            });
        }
        setColumns(cols);
        setDataset(rows);
    }, []);

    return (
        <div className="relative w-[100%] h-[100%]">
            <SlickgridReact
                gridId={"react-table-demo"}
                columns={columns}
                dataset={dataset}
                options={gridOptions}
                onReactGridCreated={onGridCreated}
            />
        </div>
    );
});

export default TableChartReactComponent;


