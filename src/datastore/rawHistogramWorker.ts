
type HistogramConfig = {
    data: SharedArrayBuffer;
    min: number;
    max: number;
    bins: number;
    isInt32: boolean;
    // todo pass in related to ~background_filter
    // as it may pertain to the view, or the chart, or in general some node in a graph...
    // filteredIndices?: SharedArrayBuffer; 
}

function computeBinsCPU(props: HistogramConfig) {
    const { isInt32, data, min, max, bins } = props;
    const arrType = isInt32 ? Int32Array : Float32Array;
    const dataArray = new arrType(data);
    const hist = new Array(bins).fill(0);
    const binWidth = (max - min) / bins;
    for (let i = 0; i < dataArray.length; i++) {
        const bin = Math.floor((dataArray[i] - min) / binWidth);
        if (bin >= 0 && bin < bins) {
            hist[bin]++;
        }
    }
    return hist;
}

// WebGPU compute shader code - with help from ChatGPT. (although it's very naive at the moment)
// thinking about a version adapting to other Dimensions...
// we'd have some boilerplate code that would be the same for all dimensions, arranging data within workgroups etc
// but the predicate logic would be different for each type of dimension
// write a bit of wgsl predicate logic to determine whether a given array element is filtered
// for each type of dimension, that would be inserted into the boilerplate shader code
// meaning that as we refine the data sharing within workgroups etc, it should hopefully
// not require changing the predicate logic, just the data handling logic



// !inject `struct Histogram` with BINS into the shader code when creating the shader module
const shaderCode = /* wgsl */ `
// struct Histogram {
//     counts: array<atomic<u32>, BINS> // Declare the array with a fixed size
// };
// struct InputData {
//     data: array<f32> // declare the array type based on isInt32
// };
struct Params {
    min: f32,
    binWidth: f32,
    bins: u32
};



@group(0) @binding(0) var<uniform> params : Params;
@group(0) @binding(1) var<storage, read> inputData : InputData;
@group(0) @binding(2) var<storage, read_write> histogram : Histogram;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) global_id : vec3<u32>) {
    let index = global_id.x;

    // Guard against out-of-bounds access to input data array
    if (index >= arrayLength(&inputData.data)) {
        return;
    }

    let value = inputData.data[index];
    let bin = u32(floor((value - params.min) / params.binWidth));

    // Ensure the computed bin is within range
    if (bin < params.bins) {
        // Perform atomic addition to ensure thread safety
        atomicAdd(&histogram.counts[bin], 1u);
    }
}
`;

async function runComputeShader(
    device: GPUDevice,
    dataBuffer: GPUBuffer,
    histBuffer: GPUBuffer,
    min: number,
    max: number,
    bins: number,
    dataLength: number,
    isInt32: boolean
): Promise<void> {
    // Create shader module
    const shaderModule = device.createShaderModule({
        code: `
        struct Histogram {
            counts: array<atomic<u32>, ${bins}> // Declare the array with a fixed size
        };
        struct InputData {
            data: array<${isInt32 ? 'i32' : 'f32'}> // declare the array type based on isInt32
        };

        ${shaderCode}
        `,
        label: 'histogram-compute-shader',
    });

    // Create bind group layout and pipeline layout
    const bindGroupLayout = device.createBindGroupLayout({
        entries: [
            { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
            { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
        ],
        label: 'histogram-bind-group-layout',
    });

    const pipelineLayout = device.createPipelineLayout({
        bindGroupLayouts: [bindGroupLayout],
        label: 'histogram-pipeline-layout',
    });

    const maxWorkgroupsPerDimension = device.limits.maxComputeWorkgroupSizeX;
    const workgroupSize = 64;

    // Create compute pipeline
    const computePipeline = device.createComputePipeline({
        layout: pipelineLayout,
        compute: {
            module: shaderModule,
            entryPoint: 'main',
        },
        label: 'histogram-compute-pipeline',
    });

    // Create the uniform buffer for the parameters
    // todo: think about a nicer abstraction for this
    const paramsBuffer = device.createBuffer({
        size: 12, // 3 * 4 bytes (min, binWidth, bins)
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
        label: 'histogram-params-buffer',
    });

    // Calculate the bin width and write parameters
    const binWidth = (max - min) / bins;
    const params = new Float32Array([min, binWidth, bins]);
    device.queue.writeBuffer(paramsBuffer, 0, params.buffer);

    // Create bind group
    const bindGroup = device.createBindGroup({
        layout: bindGroupLayout,
        entries: [
            { binding: 0, resource: { buffer: paramsBuffer } },
            { binding: 1, resource: { buffer: dataBuffer } },
            { binding: 2, resource: { buffer: histBuffer } },
        ],
        label: 'histogram-bind-group',
    });

    // Dispatch the compute shader
    const commandEncoder = device.createCommandEncoder({ label: 'histogram-command-encoder' });
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


async function computeBinsGPU(props: HistogramConfig) {
    const { isInt32, data, min, max, bins } = props;

    // Initialize WebGPU
    if (!navigator.gpu) {
        throw new Error('WebGPU is not supported by your browser.');
    }

    const adapter = await navigator.gpu.requestAdapter();
    const device = await adapter.requestDevice({label: 'histogram-device'});

    // Create GPU buffer for input data
    const arrType = isInt32 ? Int32Array : Float32Array;
    const dataArray = new arrType(data);

    // Create a buffer for the data
    const dataBuffer = device.createBuffer({
        size: dataArray.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
        label: 'data-buffer',
    });

    device.queue.writeBuffer(dataBuffer, 0, dataArray);

    // Create a buffer for the histogram (output)
    const histBuffer = device.createBuffer({
        size: bins * 4, // 4 bytes per bin (assuming 32-bit integers)
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
        label: 'histogram-buffer',
    });

    // Create compute shader (next step)
    // Set up compute pipeline (next step) 
    const len = isInt32 ? dataArray.length / 4 : dataArray.length;
    // Get the maximum buffer binding size from the device's limits
    const maxBindingSize = device.limits.maxStorageBufferBindingSize / 2048; //!! 2048 is a magic number
    const elementSize = isInt32 ? 4 : 4; // 4 bytes per element for both Int32Array and Float32Array
    const maxElementsPerChunk = Math.floor(maxBindingSize / elementSize);

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
        await runComputeShader(
            device, dataChunkBuffer, histBuffer, min, max, bins, chunkLength, isInt32
        );
    }

    // Read results from GPU and return
    const readBuffer = device.createBuffer({
        size: bins * 4, // to store the result from GPU
        usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ,
        label: 'histogram-read-buffer',
    });

    // Copy results from GPU to CPU-readable buffer
    const commandEncoder = device.createCommandEncoder({label: 'histogram-read-encoder'});
    commandEncoder.copyBufferToBuffer(histBuffer, 0, readBuffer, 0, bins * 4);
    const commands = commandEncoder.finish();
    device.queue.submit([commands]);

    // Map the result buffer to read
    await readBuffer.mapAsync(GPUMapMode.READ);
    const resultArray = new Uint32Array(readBuffer.getMappedRange());
    return Array.from(resultArray);
};


self.onmessage = async (event: MessageEvent<HistogramConfig>) => {
    // this logic is ok for MDV columns as of this writing
    // but what about viv raster data for example?
    // TODO GPU chunking of data for large arrays
    
    // todo: add options for CPU/GPU mode etc if needed
    // const data = await computeBinsGPU(event.data);
    const [histCpu, histGpu] = await Promise.all([
        // data, data
        computeBinsCPU(event.data),
        computeBinsGPU(event.data),
    ]);
    self.postMessage({histCpu, histGpu});
    // const hist = await computeBinsGPU(event.data);
    // self.postMessage(hist);
}
