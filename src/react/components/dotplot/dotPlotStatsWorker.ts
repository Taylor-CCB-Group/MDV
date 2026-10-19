import { type DotPlotStatsRequest, computeDotPlotStats } from "./dotPlotStats";

self.onmessage = (event: MessageEvent<DotPlotStatsRequest>) => {
    self.postMessage(computeDotPlotStats(event.data));
};
