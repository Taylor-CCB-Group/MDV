import type { DataType, LoadedDataColumn } from "@/charts/charts";
import type DataStore from "@/datastore/DataStore";

class SlickGridDataProvider {
  private columns: LoadedDataColumn<DataType>[];
  private indices: Uint32Array;
  private dataStore: DataStore;
  private includeIndex: boolean;
  constructor(dataStore: DataStore, columns: LoadedDataColumn<DataType>[], indices: Uint32Array, includeIndex = false) {
    this.columns = columns;
    this.indices = indices;
    this.dataStore = dataStore;
    this.includeIndex = includeIndex;
  }

  getLength(): number {
    return this.indices.length;
  }

  getItem(index: number): any {
    if (index < 0 || index >= this.indices.length) {
      return null;
    }

    const i = this.indices[index]
    const item: any = {};

    if (this.includeIndex) {
      item.__index__ = i + 1;
    }

    for (const col of this.columns) {
      //! Make use col.getValue instead of dataStore directly here to get the value. Note that getRowAsObject might not be optimal for large data as stated in it's description.
      item[col.name] = this.dataStore.getRowText(i, col.field);
    }

    return item;
  }

  getItemMetadata(_index: number) {
    return null;
  }

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