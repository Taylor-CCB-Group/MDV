import type { ViewState as SpatialCanvasViewState } from "@spatialdata/vis";
import type { OrbitViewState, OrthographicViewState } from "deck.gl";

export function toSpatialViewState(
    viewState: OrthographicViewState | OrbitViewState | null | undefined,
): SpatialCanvasViewState | null {
    if (!viewState || !Array.isArray(viewState.target)) return null;
    const [x, y] = viewState.target;
    if (typeof x !== "number" || typeof y !== "number") return null;
    const zoom = typeof viewState.zoom === "number" ? viewState.zoom : 0;
    return {
        target: [x, y],
        zoom,
    };
}

export function toMdvViewState(
    viewState: SpatialCanvasViewState,
    previousViewState: OrthographicViewState | OrbitViewState | null | undefined,
): OrthographicViewState {
    const [x, y] = viewState.target;
    const previousTarget = Array.isArray(previousViewState?.target)
        ? previousViewState.target
        : [0, 0, 0];
    return {
        ...(previousViewState ?? {}),
        target: [x, y, previousTarget[2] ?? 0],
        zoom: viewState.zoom,
    };
}
