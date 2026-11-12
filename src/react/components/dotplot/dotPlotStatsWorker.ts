import { computeDotPlotStats, type DotPlotStatsRequest } from "./dotPlotStats";

self.onmessage = (event: MessageEvent<DotPlotStatsRequest>) => {
    self.postMessage(computeDotPlotStats(event.data));
};
