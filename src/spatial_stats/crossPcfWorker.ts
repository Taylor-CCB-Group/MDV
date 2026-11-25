import { type CrossPcfRequest, crossPcfInDomain } from "./crossPcfCompute";

self.onmessage = (e: MessageEvent<CrossPcfRequest>) => {
    try {
        const result = crossPcfInDomain(e.data);
        self.postMessage({ result }, { transfer: [result.g.buffer, result.pairs.buffer] });
    } catch (error) {
        self.postMessage({ error: error instanceof Error ? error.message : String(error) });
    }
};
