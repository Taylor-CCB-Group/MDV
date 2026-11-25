import { computeMultiLineDensities, type MultiLineDensityRequest } from "./multiLineStats";

self.onmessage = (event: MessageEvent<MultiLineDensityRequest>) => {
    self.postMessage(computeMultiLineDensities(event.data));
};
