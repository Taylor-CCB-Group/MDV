
type HistogramConfig = {
    data: SharedArrayBuffer;
    min: number;
    max: number;
    bins: number;
    arrayType: "float32" | "int32" | "uint32";
    byteOffset: number;
    length: number;
    // todo pass in related to ~background_filter
    // as it may pertain to the view, or the chart, or in general some node in a graph...
    // filteredIndices?: SharedArrayBuffer; 
}

self.onmessage = async (event: MessageEvent<HistogramConfig>) => {
    const { arrayType, data, min, max, bins, byteOffset, length } = event.data;
    // this logic is ok for MDV columns as of this writing
    // but what about viv raster data for example?
    const dataArray =
        arrayType === "int32"
            ? new Int32Array(data, byteOffset, length)
            : arrayType === "uint32"
              ? new Uint32Array(data, byteOffset, length)
              : new Float32Array(data, byteOffset, length);
    const hist = new Array(bins).fill(0);
    const binWidth = (max - min) / bins;
    for (let i = 0; i < dataArray.length; i++) {
        const value = dataArray[i];
        const bin = value === max ? bins - 1 : Math.floor((value - min) / binWidth);
        if (bin >= 0 && bin < bins) {
            hist[bin]++;
        }
    }
    self.postMessage(hist);
}
