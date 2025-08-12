import type { DataColumn, DataSource } from "@/charts/charts";
import type { DataLoader } from "./DataLoaderUtil";


export default function getZarrDataLoader(datasources: DataSource[], root: string): DataLoader {
  throw new Error("Not implemented");
  // return {
    // function: async (columns: DataColumn<any>[], dataSource: string, size: number) => {
    //     const zarr = datasources.find(ds => ds.name === dataSource)?.zarr; 
    // },
    // viewLoader: async (view: string) => views[view],
    // rowDataLoader: loadRowDataStatic,
    // binaryDataLoader: loadBinaryDataStatic,
  // }
};