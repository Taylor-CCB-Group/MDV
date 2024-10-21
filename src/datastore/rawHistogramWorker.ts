
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

const shaderCode = /* wgsl */ `
struct Params {
    min: f32,
    binWidth: f32,
    bins: u32
};

struct InputData {
    data: array<f32>
};

struct Histogram {
    counts: array<atomic<u32>, 1> // Declare the array with a fixed size
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
    dataLength: number
): Promise<void> {
    // Create shader module
    const shaderModule = device.createShaderModule({
        code: shaderCode,
    });

    // Create bind group layout and pipeline layout
    const bindGroupLayout = device.createBindGroupLayout({
        entries: [
            { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
            { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
        ],
    });

    const pipelineLayout = device.createPipelineLayout({
        bindGroupLayouts: [bindGroupLayout],
    });

    // Create compute pipeline
    const computePipeline = device.createComputePipeline({
        layout: pipelineLayout,
        compute: {
            module: shaderModule,
            entryPoint: 'main',
        },
    });

    // Create the uniform buffer for the parameters
    const paramsBuffer = device.createBuffer({
        size: 12, // 3 * 4 bytes (min, binWidth, bins)
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
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
    });

    // Dispatch the compute shader
    const commandEncoder = device.createCommandEncoder();
    const passEncoder = commandEncoder.beginComputePass();
    passEncoder.setPipeline(computePipeline);
    passEncoder.setBindGroup(0, bindGroup);
    passEncoder.dispatchWorkgroups(Math.ceil(dataLength / 64)); // Now we pass `dataLength` for dispatch
    passEncoder.end();

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
    const device = await adapter.requestDevice();

    // Create GPU buffer for input data
    const arrType = isInt32 ? Int32Array : Float32Array;
    const dataArray = new arrType(data);

    // Create a buffer for the data
    const dataBuffer = device.createBuffer({
        size: dataArray.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });

    device.queue.writeBuffer(dataBuffer, 0, dataArray);

    // Create a buffer for the histogram (output)
    const histBuffer = device.createBuffer({
        size: bins * 4, // 4 bytes per bin (assuming 32-bit integers)
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
    });

    // Create compute shader (next step)
    // Set up compute pipeline (next step) 
    await runComputeShader(device, dataBuffer, histBuffer, min, max, bins, dataArray.length / 4);

    // Read results from GPU and return
    const readBuffer = device.createBuffer({
        size: bins * 4, // to store the result from GPU
        usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ,
    });

    // Copy results from GPU to CPU-readable buffer
    const commandEncoder = device.createCommandEncoder();
    commandEncoder.copyBufferToBuffer(histBuffer, 0, readBuffer, 0, bins * 4);
    const commands = commandEncoder.finish();
    device.queue.submit([commands]);

    // Map the result buffer to read
    await readBuffer.mapAsync(GPUMapMode.READ);
    const resultArray = new Int32Array(readBuffer.getMappedRange());
    return Array.from(resultArray);
};


self.onmessage = async (event: MessageEvent<HistogramConfig>) => {
    const { isInt32, data, min, max, bins } = event.data;
    // this logic is ok for MDV columns as of this writing
    // but what about viv raster data for example?
    const hist = await computeBinsGPU({ isInt32, data, min, max, bins });
    self.postMessage(hist);
}
