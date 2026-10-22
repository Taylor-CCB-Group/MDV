import { type FilterRequest, handleFilterRequest } from "./rangeFilterEval";

// Evaluates a RangeDimension filter spec into a shared include mask; see rangeFilterEval.ts.
self.onmessage = (event: MessageEvent<FilterRequest>) => {
    self.postMessage(handleFilterRequest(event.data));
};
