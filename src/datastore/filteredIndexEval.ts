/**
 * Lists the rows that pass a DataStore's filter. Shared by `filteredIndexWorker.ts` and the main thread, which
 * does the same work itself if the worker can't be used.
 */

export type FilteredIndexRequest = {
    /** echoed in the response: the worker is shared, so replies have to be matched to requests */
    id: number;
    /** the store's per-row filter state: 0 means the row passes all filters */
    filterBuffer: SharedArrayBuffer;
    /** a `Uint32Array` buffer with room for every row that passed when the request was made */
    outputBuffer: SharedArrayBuffer;
};

export type FilteredIndexResponse = { type: "done"; id: number };

/** Writes the index of every row that passes the filter into `output`, in ascending order. */
export function writeFilteredIndices(filterArray: Uint8Array, output: Uint32Array) {
    for (let i = 0, j = 0; i < filterArray.length; i++) {
        if (filterArray[i] === 0) {
            output[j++] = i;
        }
    }
}

/** The worker's message handler: fills the request's output buffer and returns the reply to post. */
export function handleFilteredIndexRequest(request: FilteredIndexRequest): FilteredIndexResponse {
    //todo atomics.
    writeFilteredIndices(new Uint8Array(request.filterBuffer), new Uint32Array(request.outputBuffer));
    return { type: "done", id: request.id };
}
