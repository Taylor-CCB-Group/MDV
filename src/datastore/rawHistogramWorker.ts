
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

self.onmessage = async (event: MessageEvent<HistogramConfig>) => {
    const { isInt32, data, min, max, bins } = event.data;
    // this logic is ok for MDV columns as of this writing
    // but what about viv raster data for example?
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
    self.postMessage(hist);
}