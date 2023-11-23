
onmessage = function (e: MessageEvent<{ filterBuffer: SharedArrayBuffer, outputBuffer: SharedArrayBuffer }>) {
    //todo atomics.
    const filterArray = new Uint8Array(e.data.filterBuffer);
    const outputArray = new Uint32Array(e.data.outputBuffer);
    for (let i=0, j=0; i<filterArray.length; i++) {
        if (filterArray[i] === 0) {
            outputArray[j++] = i;
        }
    }
    postMessage({ type: "done" });
};