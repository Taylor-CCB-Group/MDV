// biome-ignore lint/suspicious/noGlobalAssign: <explanation>
onmessage = (
    e: MessageEvent<{
        filterBuffer: SharedArrayBuffer;
        outputBuffer: SharedArrayBuffer;
        id: number;
    }>,
) => {
    //todo atomics.
    const filterArray = new Uint8Array(e.data.filterBuffer);
    const outputArray = new Uint32Array(e.data.outputBuffer);
    for (let i = 0, j = 0; i < filterArray.length; i++) {
        if (filterArray[i] === 0) {
            outputArray[j++] = i;
        }
    }
    postMessage({ type: "done", id: e.data.id });
};
