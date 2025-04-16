import type { TypedArray } from "d3";

export type HistogramParams = {
    min: number;
    max: number;
    bins: number;
    // we should probably have more than two supported data types...
    isInt32: boolean;
}

export type HistogramInput = {
    data: TypedArray;
} & HistogramParams;

export type HistogramMessage = {
    data: SharedArrayBuffer;
    // todo pass in related to ~background_filter
    // as it may pertain to the view, or the chart, or in general some node in a graph...
    // filteredIndices?: SharedArrayBuffer;
} & HistogramParams;

self.onmessage = async (event: MessageEvent<HistogramMessage>) => {
    const { isInt32, data, min, max, bins } = event.data;
    // this logic is ok for MDV columns as of this writing
    // but what about viv raster data for example?
    const arrType = isInt32 ? Int32Array : Float32Array;
    //@ts-ignore !not sure why there's a discrepancy here between tsc & language server in vscode
    const dataArray = new arrType(data);
    const hist = new Array(bins).fill(0);
    // todo - consider using d3 scale for bins
    const binWidth = (max - min) / bins;
    for (let i = 0; i < dataArray.length; i++) {
        const bin = Math.floor((dataArray[i] - min) / binWidth);
        if (bin >= 0 && bin < bins) {
            hist[bin]++;
        }
    }
    self.postMessage(hist);
}
