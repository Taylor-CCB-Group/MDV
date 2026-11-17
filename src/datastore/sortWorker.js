import { sortOrder } from "./sortColumns";

// biome-ignore lint/suspicious/noGlobalAssign: a worker's message handler is the global onmessage
onmessage = (e) => {
    const { orderBuffer, columns } = e.data;
    sortOrder(new Uint32Array(orderBuffer), columns);
    postMessage("done");
};
