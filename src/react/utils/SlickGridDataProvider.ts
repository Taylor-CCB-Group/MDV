import type { DataType, LoadedDataColumn } from "@/charts/charts";
import type DataStore from "@/datastore/DataStore";

class SlickGridDataProvider {
  private columns: LoadedDataColumn<DataType>[];
  private indices: Uint32Array;
  private dataStore: DataStore;
  private includeIndex: boolean;
  private selectedRowIds: Set<number>;
  constructor(dataStore: DataStore, columns: LoadedDataColumn<DataType>[], indices: Uint32Array, includeIndex = false) {
    this.columns = columns;
    this.indices = indices;
    this.dataStore = dataStore;
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

    // Use getRowAsObject exactly like old DataModel.getItem() does
    const columnFields = this.columns.map(c => c.field);
    const item: any = this.dataStore.getRowAsObject(dataIndex, columnFields);
    if (!item) return null;
    // Add id for SlickGrid
    item.id = index;

    // Override __index__ if we're showing the index column
    if (this.includeIndex) {
      item.__index__ = dataIndex + 1;
    }

    return item;
  }

  getItemMetadata(_index: number) {
    return null;
  }

  // Get the data index for a grid row (for highlighting)
  getDataIndex(gridRow: number): number {
    return this.indices[gridRow];
  }

  // // Find grid row for a data index (for external highlighting)
  // findGridRow(dataIndex: number): number {
  //   return this.indices.indexOf(dataIndex);
  // }

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
        .filter(row => row >= 0 && row < this.indices.length)
        .map(row => {
            const item = this.getItem(row);
            return item?.id ?? row;
        });
}

  // // Map row IDs to row indices
  // mapIdsToRows(ids: number[]): number[] {
  //   return ids.filter(id => id >= 0 && id < this.indices.length);
  // }

  // // Map row indices to row IDs
  // mapRowsToIds(rows: number[]): number[] {
  //   return rows.filter(row => row >= 0 && row < this.indices.length);
  // }

  // getDataIndex(gridRow: number) {
  //     return this.indices[gridRow];
  // }

  // findGridRow(dataIndex: number) {
  //     return this.indices.indexOf(dataIndex);
  // }

  // Avoid calling the default methods for sorting and filtering
  sort() {

  }

  reSort() {

  }

  refresh() {

  }
}

export default SlickGridDataProvider;