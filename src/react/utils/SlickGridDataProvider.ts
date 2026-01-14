import type { DataType, LoadedDataColumn } from "@/charts/charts";

class SlickGridDataProvider {
    private columns: LoadedDataColumn<DataType>[];
    private indices: Uint32Array;
    private includeIndex: boolean;
    private selectedRowIds: Set<number>;
    constructor(columns: LoadedDataColumn<DataType>[], indices: Uint32Array, includeIndex = false) {
        this.columns = columns;
        this.indices = indices;
        this.includeIndex = includeIndex;
        this.selectedRowIds = new Set();
    }

    getLength(): number {
        return this.indices.length;
    }

    getItem(index: number): any {
        if (index < 0 || index >= this.indices.length) {
            return null;
        }

        const dataIndex = this.indices[index];

        const item: any = {};
        for (const column of this.columns) {
            if (column) item[column.field] = column.getValue(dataIndex);
        }

        if (!item) return null;
        // Add id for SlickGrid
        item.id = dataIndex;

        // Override __index__ if we're showing the index column
        if (this.includeIndex) {
            item.__index__ = dataIndex + 1;
        }

        return item;
    }

    getItemMetadata(_index: number) {
        return null;
    }

    // Required by ResizerService for resize by cell content
    // Return all the items of data provider
    getItems(): any[] {
        const items: any[] = [];
        for (let i = 0; i < this.indices.length; i++) {
            const item = this.getItem(i);
            if (item) items.push(item);
        }

        return items;
    }

    // Get the data index for a grid row (for highlighting)
    getDataIndex(gridRow: number): number {
        if (gridRow < 0 || gridRow >= this.indices.length) {
            return -1;
        }
        return this.indices[gridRow] ?? -1;
    }

    // Required by GridStateService for row selection
    getAllSelectedIds(): number[] {
        return Array.from(this.selectedRowIds);
    }

    // Required by GridStateService for filtered selection
    getAllSelectedFilteredIds(): number[] {
        return Array.from(this.selectedRowIds);
    }

    // Set selected row IDs
    setSelectedIds(ids: number[]) {
        this.selectedRowIds = new Set(ids);
    }

    mapRowsToIds(rows: number[]): number[] {
        return rows
            .filter((row) => row >= 0 && row < this.indices.length)
            .map((row) => {
                // Return stable data index (not grid row position)
                return this.indices[row];
            });
    }

    // Avoid calling the default methods for sorting and filtering
    sort() {}

    reSort() {}

    refresh() {}
}

export default SlickGridDataProvider;
