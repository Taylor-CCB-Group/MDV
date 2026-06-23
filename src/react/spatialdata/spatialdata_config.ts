export type SpatialDataSerializableViewState = {
    target: number[];
    zoom?: number;
};

export function removeSpatialDataRootViv<T extends object>(config: T): T {
    delete (config as T & { viv?: unknown }).viv;
    return config;
}

export function toSerializableSpatialDataViewState(
    viewState: SpatialDataSerializableViewState | null | undefined,
): SpatialDataSerializableViewState | null {
    if (!viewState?.target) return null;
    return {
        target: viewState.target,
        zoom: viewState.zoom,
    };
}
