import { DataColumn, DataType } from "@/charts/charts";
import DataStore from "./DataStore";

// 
export type ColumnGPU<T extends DataType> = {
    column: DataColumn<T>;
    gpuBuffer: GPUBuffer;
    device: GPUDevice;
}

export async function getColumnBuffer(column: DataColumn<DataType>) {
    // if it is indeed basically a singleton, we shouldn't need to think about this DataStoreGPU business
    // return something like the above?
    return (await DataStoreGPU.getDataStoreGPU()).getColumnBuffer(column);
}

export default class DataStoreGPU { //extends DataStore?
    device: GPUDevice;
    columnBuffers: Map<DataColumn<DataType>, GPUBuffer> = new Map();
    constructor(device: GPUDevice) {
        // actually, the expectation that there should be a one-to-one mapping between
        // DataStore and GPUDevice implies that we don't do GPU operations on multiple DataStores
        // so that doesn't seem like the right design.
        this.device = device;
    }
    async getColumnBuffer(column: DataColumn<DataType>) {
        if (!column.data) {
            // await loadColumn(column);
            throw new Error("Column data is not available (we should load it here)");
        }
        const { device } = this;
        const { data } = column;
        // todo data type related 
        const dataBuffer = device.createBuffer({
            size: data.byteLength,
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
        });
        device.queue.writeBuffer(dataBuffer, 0, data);
        // we probably want to return our own type here, which includes the `dataBuffer`,
        // but also some more context/metadata about the column
        // Is there any reason columns shouldn't have a reference to DataStore?
        return dataBuffer;
    }
    // not sure how keen I am on this.
    //! we probably could avoid leaking resources all over the place if we keep things more local
    // worth thinking about lifetimes & ownership, releasing resources properly etc.
    static singleton?: DataStoreGPU
    static async getDataStoreGPU() {
        if (DataStoreGPU.singleton) return DataStoreGPU.singleton;
        if (!navigator.gpu) {
            throw new Error("WebGPU not supported in this browser.");
        }
        const adapter = await navigator.gpu.requestAdapter();
        const device = await adapter?.requestDevice({ label: "datastore-gpu" });
        if (!device) {
            throw new Error("No device available");
        }
        DataStoreGPU.singleton = new DataStoreGPU(device);
        return DataStoreGPU.singleton;
    }
}