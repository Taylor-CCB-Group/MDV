export type SpatialDataSerializableViewState = {
    target: number[];
    zoom?: number;
};

export function toSerializableSpatialDataViewState(
    viewState: SpatialDataSerializableViewState | null | undefined,
): SpatialDataSerializableViewState | null {
    if (!viewState?.target) return null;
    return {
        target: viewState.target,
        zoom: viewState.zoom,
    };
}
