import type { LoadedDataColumn, NumberDataType } from "@/charts/charts";
import Dimension from "./Dimension";
import DataStoreGPU from "./DataStoreGPU";
import DataStore from "./DataStore";

export type RangeArgs = {
    min: number,
    max: number,
};

const WG_SIZE = 64;

// TODO parameterise data type:
// for multi-parameter filters, on a per-parameter basis...
// so: general way to represent columns in WebGPU for other purposes
// design will be clearer when we have a clearer idea of how they will be used/bound etc,.
const rangeShaderDeclCode = /* wgsl */`
struct Params {
    min: f32,
    max: f32
};
//this type will depend on the type of column data,
//and also potentially multiple columns...
struct InputData {
    data: array<f32>
};
@group(0) @binding(2) var<storage, read> inputData : InputData;
//@group(0) @binding(0) var<uniform> params : Params; //probably common rather than injected
// do we pass an index or a value? do we care? probably can't use before declaring.
fn predicate(index: u32) -> bool {
    let value = &inputData.data[index];
    return value >= params.min && value <= params.max;
}
`;

const rangeShaderCode = /* wgsl */`
${rangeShaderDeclCode}

struct FilterState {
    localFilter: array<u8>, //RW
    backgroundFilter: array<u8>, //R
    filterOutput: array<u8>, //
}

@group(0) @binding(0) var<uniform> params : Params;
// what we actually need here is read-write FilterState, not sure about semantics
//@group(0) @binding(1) var<uniform> filterState : FilterState; //?
@group(0) @binding(1) var<storage, write> filterOutput : array<u8>;

@compute @workgroup_size(${WG_SIZE})
fn main(@builtin(global_invocation_id) global_id : vec3<u32>) {
    let i = global_id.x;
    if (i >= arrayLength(&inputData.data)) {
        return;
    }
    // this should mimic the behaviour of the CPU version
    // filterSize - atomic, or on CPU after?
    let value = !predicate(i);
    // JS for ref
    // if (value) {
    //     if (localFilter[i] === 0) {
    //         if (++filter[i] === 1) {
    //             parent.filterSize--;
    //         }
    //     }
    //     localFilter[i] = 1;
    // } else {
    //     if (localFilter[i] === 1) {
    //         if (--filter[i] === 0) {
    //             parent.filterSize++;
    //         }
    //     }
    //     localFilter[i] = 0;
    // }

    &filterOutput[i] = u8(value);
}
/**
 * Copilot generated pseudo code for multiple passes...
 */

// First pass: Compute contributions
@compute @workgroup_size(64)
fn compute_contributions(@builtin(global_invocation_id) global_id : vec3<u32>) {
    let index = global_id.x;
    if (index >= data_length) {
        return;
    }
    let contribution = compute_filter_contribution(data[index]);
    contributions[index] = contribution;
}

// Second pass: Parallel reduction
@compute @workgroup_size(64)
fn reduce_contributions(@builtin(global_invocation_id) global_id : vec3<u32>) {
    let local_id = global_id.x % 64;
    let group_id = global_id.x / 64;

    // Load contributions into shared memory
    shared local_contributions: array<f32, 64>;
    local_contributions[local_id] = contributions[global_id.x];
    workgroupBarrier();

    // Perform reduction within the workgroup
    var stride = 32u;
    while (stride > 0u) {
        if (local_id < stride) {
            local_contributions[local_id] += local_contributions[local_id + stride];
        }
        stride = stride / 2u;
        workgroupBarrier();
    }

    // Write the result of this workgroup to the output buffer
    if (local_id == 0u) {
        reduced[group_id] = local_contributions[0];
    }
}

// Final reduction pass (if needed)
@compute @workgroup_size(64)
fn final_reduce(@builtin(global_invocation_id) global_id : vec3<u32>) {
    // Similar reduction logic to sum up the reduced buffer into a single value
}

// Third pass: Apply filter
@compute @workgroup_size(64)
fn apply_filter(@builtin(global_invocation_id) global_id : vec3<u32>) {
    let index = global_id.x;
    if (index >= data_length) {
        return;
    }
    let result = apply_filter_logic(data[index], filterSize);
    output[index] = result;
}
`;

/** copilot version */
async function runWebGPUFilter(data: Float32Array) {
    // Step 1: Initialize WebGPU
    const adapter = await navigator.gpu.requestAdapter();
    const device = await adapter?.requestDevice();
    if (!device) throw new Error("Failed to initialize WebGPU");

    const queue = device.queue;

    // Create buffers
    const dataBuffer = device.createBuffer({
        size: data.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
        mappedAtCreation: true,
    });
    new Float32Array(dataBuffer.getMappedRange()).set(data);
    dataBuffer.unmap();

    const contributionsBuffer = device.createBuffer({
        size: data.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST,
    });

    const reducedBuffer = device.createBuffer({
        size: 256 * Float32Array.BYTES_PER_ELEMENT, // Assuming max 256 workgroups
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST,
    });

    const filterSizeBuffer = device.createBuffer({
        size: Float32Array.BYTES_PER_ELEMENT,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST,
    });

    const outputBuffer = device.createBuffer({
        size: data.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
    });

    const resultBuffer = device.createBuffer({
        size: data.byteLength,
        usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ,
    });

    // Step 2: Load and compile shaders
    const shaderModule = device.createShaderModule({
        code: `
        // WGSL shader code from the previous example goes here
        `,
    });

    // Step 3: Create compute pipelines
    const computeContributionsPipeline = device.createComputePipeline({
        layout: "auto",
        compute: { module: shaderModule, entryPoint: "compute_contributions" },
    });

    const reduceContributionsPipeline = device.createComputePipeline({
        layout: "auto",
        compute: { module: shaderModule, entryPoint: "reduce_contributions" },
    });

    const applyFilterPipeline = device.createComputePipeline({
        layout: "auto",
        compute: { module: shaderModule, entryPoint: "apply_filter" },
    });

    // Step 4: Create bind groups
    const computeContributionsBindGroup = device.createBindGroup({
        layout: computeContributionsPipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: dataBuffer } },
            { binding: 1, resource: { buffer: contributionsBuffer } },
        ],
    });

    const reduceContributionsBindGroup = device.createBindGroup({
        layout: reduceContributionsPipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: contributionsBuffer } },
            { binding: 1, resource: { buffer: reducedBuffer } },
        ],
    });

    const applyFilterBindGroup = device.createBindGroup({
        layout: applyFilterPipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: dataBuffer } },
            { binding: 1, resource: { buffer: filterSizeBuffer } },
            { binding: 2, resource: { buffer: outputBuffer } },
        ],
    });

    // Step 5: Create command encoder and dispatch passes
    const commandEncoder = device.createCommandEncoder();

    // First pass: Compute contributions
    {
        const passEncoder = commandEncoder.beginComputePass();
        passEncoder.setPipeline(computeContributionsPipeline);
        passEncoder.setBindGroup(0, computeContributionsBindGroup);
        passEncoder.dispatchWorkgroups(Math.ceil(data.length / 64));
        passEncoder.end();
    }

    // Second pass: Reduce contributions
    {
        const passEncoder = commandEncoder.beginComputePass();
        passEncoder.setPipeline(reduceContributionsPipeline);
        passEncoder.setBindGroup(0, reduceContributionsBindGroup);
        passEncoder.dispatchWorkgroups(256); // Assuming max 256 workgroups
        passEncoder.end();
    }

    // Third pass: Apply filter
    {
        const passEncoder = commandEncoder.beginComputePass();
        passEncoder.setPipeline(applyFilterPipeline);
        passEncoder.setBindGroup(0, applyFilterBindGroup);
        passEncoder.dispatchWorkgroups(Math.ceil(data.length / 64));
        passEncoder.end();
    }

    // Copy output buffer to result buffer
    commandEncoder.copyBufferToBuffer(outputBuffer, 0, resultBuffer, 0, data.byteLength);

    // Submit commands
    queue.submit([commandEncoder.finish()]);

    // Step 6: Read back results
    await resultBuffer.mapAsync(GPUMapMode.READ);
    const resultArray = new Float32Array(resultBuffer.getMappedRange());
    console.log("Filtered Data:", resultArray);
    resultBuffer.unmap();
}



async function runComputeShader(
    device: GPUDevice,
    dataBuffer: GPUBuffer,
    filterBuffer: GPUBuffer,
    min: number,
    max: number,
    dataLength: number,
    isInt32: boolean
) {
    const shaderModule = device.createShaderModule({
        code: rangeShaderCode
    });

    const maxWorkgroupsPerDimension = device.limits.maxComputeWorkgroupSizeX;
    const workgroupSize = 64;

    // Create compute pipeline
    const computePipeline = device.createComputePipeline({
        layout: 'auto',
        compute: {
            module: shaderModule,
            entryPoint: 'main',
        },
    });

    // Create the uniform buffer for the parameters
    // todo: think about a nicer abstraction for this
    const paramsBuffer = device.createBuffer({
        size: 8, // 2 * 4 bytes (min, max)
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    // Calculate the bin width and write parameters
    const params = new Float32Array([min, max]);
    device.queue.writeBuffer(paramsBuffer, 0, params.buffer);

    // Create bind group
    const bindGroup = device.createBindGroup({
        layout: computePipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: paramsBuffer } },
            { binding: 1, resource: { buffer: filterBuffer } },
            { binding: 2, resource: { buffer: dataBuffer } },
        ],
    });

    // Dispatch the compute shader
    const commandEncoder = device.createCommandEncoder();
    // Dispatch in chunks if needed
    // !this is now giving wrong results for int32 data?
    /// (also it's actually way slower than the CPU version - we need a non-naive shader)
    let remainingElements = dataLength;
    let offset = 0;

    while (remainingElements > 0) {
        const workgroupsForThisPass = Math.min(
            Math.ceil(remainingElements / workgroupSize),
            maxWorkgroupsPerDimension
        );

        // Begin a compute pass
        const passEncoder = commandEncoder.beginComputePass({ label: 'histogram-compute-pass' });
        passEncoder.setPipeline(computePipeline);
        passEncoder.setBindGroup(0, bindGroup);

        // Dispatch only as many workgroups as allowed
        passEncoder.dispatchWorkgroups(workgroupsForThisPass);
        passEncoder.end();

        // Update remaining elements and offset
        remainingElements -= workgroupsForThisPass * workgroupSize;
        offset += workgroupsForThisPass * workgroupSize;
    }

    // Submit commands to GPU
    device.queue.submit([commandEncoder.finish()]);
}



export default class DimensionGPU extends Dimension {
    filterBufferGpu?: GPUBuffer;
    async getGpuFilterBuffer(device: GPUDevice) {
        if (this.filterBufferGpu) return this.filterBufferGpu;
        this.filterBufferGpu = device.createBuffer({
            size: this.filterArray.byteLength,
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
        });
        return this.filterBufferGpu;
    }
    //original return type is `void`... all about side-effect of writing filterArray
    // async filterRange(range: RangeArgs, column: LoadedDataColumn<NumberDataType>) {
    async filterRange(range: RangeArgs, columnName: string) {
        // const { min, max } = range;
        // const arr = this.parent.columnIndex[columnName].data;
        // // performance seems similar to non-predicate version
        // const predicate = (i: number) => {
        //     const v = arr[i];
        //     return v >= min && v <= max && !Number.isNaN(v);
        // };
        // return this.filterPredicate({ predicate });

        //! cancel old job if running. useQuery.
        // we should just have a column object. not a string which may or may not be in columnIndex...
        const column = this.parent.columnIndex[columnName] as LoadedDataColumn<NumberDataType>;
        const dsGpu = await DataStoreGPU.getDataStoreGPU();
        const { device } = dsGpu;
        const { data } = column;
        const isInt32 = false;
        // this could probably just be a method on the column...
        // while it's experimental, we may not want to change column implementation
        
        // unused because of dataChunkBuffer being used instead?
        // const dataBuffer = dsGpu.getColumnBuffer(column);
        const filterBuffer = await this.getGpuFilterBuffer(device);

        const maxWorkgroupsPerDimension = device.limits.maxComputeWorkgroupSizeX;
        const maxBindingSize = device.limits.maxStorageBufferBindingSize / 2048;
        const elementSize = 4; //f32/i32
        const maxElementsPerChunk = Math.floor(maxBindingSize / elementSize);
        const len = data.length;

        // Calculate how many chunks we need
        const totalChunks = Math.ceil(len / maxElementsPerChunk);
        for (let chunkIndex = 0; chunkIndex < totalChunks; chunkIndex++) {
            const start = chunkIndex * maxElementsPerChunk;
            const end = Math.min(start + maxElementsPerChunk, len);

            const chunkLength = end - start;
            const dataChunkBuffer = device.createBuffer({
                size: chunkLength * elementSize,
                usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
                mappedAtCreation: true,
                label: `data-chunk-buffer-${chunkIndex}`,
            });

            // Copy the data for the chunk into the buffer
            const arrayBufferView = isInt32 ? new Int32Array(data, start * elementSize, chunkLength)
                : new Float32Array(data, start * elementSize, chunkLength);
            new Float32Array(dataChunkBuffer.getMappedRange()).set(arrayBufferView);
            dataChunkBuffer.unmap();

            // Now bind the data chunk to the shader and run the compute pass
            const { min, max } = range;
            await runComputeShader(
                device, dataChunkBuffer, filterBuffer, min, max, chunkLength, isInt32
            );
        }

        // Read results from GPU and return
        // rounding down to nearest multiple of 4 is not producing errors, or useful output
        const readSize = Math.floor(data.length / 4) * 4;
        const readBuffer = device.createBuffer({
            //! probably not the same byteLength - isn't this <u8>?
            size: readSize, // to store the result from GPU
            usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ,
        });

        // Copy results from GPU to CPU-readable buffer
        const commandEncoder = device.createCommandEncoder();
        commandEncoder.copyBufferToBuffer(filterBuffer, 0, readBuffer, 0, readSize);
        const commands = commandEncoder.finish();
        device.queue.submit([commands]);

        // Map the result buffer to read
        await readBuffer.mapAsync(GPUMapMode.READ);
        // const resultArray = new Uint8Array(readBuffer.getMappedRange());
        this.filterArray.set(new Uint8Array(readBuffer.getMappedRange()));

        
        // return Array.from(resultArray);
    }
    destroy(notify?: boolean): void {
        super.destroy(notify);
        this.filterBufferGpu?.destroy();
        this.filterBufferGpu = undefined;
    }
}
Dimension.types["range_gpu"] = DimensionGPU;