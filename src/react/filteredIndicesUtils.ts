const FILTERED_INDICES_REFRESH_EVENTS = new Set([
    "filtered",
    "data_changed",
    "data_added",
]);

export function shouldRefreshFilteredIndices(eventType: string) {
    return FILTERED_INDICES_REFRESH_EVENTS.has(eventType);
}
