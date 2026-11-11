import { type FilteredIndexRequest, handleFilteredIndexRequest } from "./filteredIndexEval";

// Lists the rows that pass a DataStore's filter into a shared buffer; see filteredIndexEval.ts.
self.onmessage = (event: MessageEvent<FilteredIndexRequest>) => {
    self.postMessage(handleFilteredIndexRequest(event.data));
};
