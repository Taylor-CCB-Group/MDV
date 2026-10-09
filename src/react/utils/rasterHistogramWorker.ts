import { computeRasterHistogram, type RasterHistogramRequest } from "./rasterHistogram";

self.onmessage = (event: MessageEvent<RasterHistogramRequest>) => {
    self.postMessage(computeRasterHistogram(event.data));
};
