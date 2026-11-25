import type { Layer, LayerProps, PickingInfo } from "@deck.gl/core";
import { getIsolatedGroupMember, IsolatedGroupLayer } from "@spatialdata/layers";
import { recordSpatialPerf, spatialPerfEnabled } from "@/react/spatialdata/perf";

/**
 * The group's id. It must not start with `density-grid-`, which would tie it to one cell
 * (see `shouldDrawInDensityGridCell`).
 */
export const SHARED_CELL_LAYERS_ID = "shared-cell-layers";

type PickingEventName = "onClick" | "onDragStart" | "onDrag" | "onDragEnd";
type PickingHandler = NonNullable<LayerProps["onClick"]>;

/**
 * deck hands a pick's events to its top-level layer, which for a layer in the group is the
 * group. Pass them on to the member the pick came through, as deck would have.
 */
function forwardToMember(name: PickingEventName): PickingHandler {
    return (info, event) => {
        const member = getIsolatedGroupMember(info);
        const handler = member?.props[name];
        return handler ? (handler.call(member, info, event) ?? false) : false;
    };
}

const forwardHover: PickingHandler = (info: PickingInfo, event) =>
    getIsolatedGroupMember(info)?.onHover(info, event) ?? false;

const forwardedEvents = {
    onHover: forwardHover,
    onClick: forwardToMember("onClick"),
    onDragStart: forwardToMember("onDragStart"),
    onDrag: forwardToMember("onDrag"),
    onDragEnd: forwardToMember("onDragEnd"),
};

/** Counts target draws, so a capture shows the group drawing once per frame (see `perf.ts`). */
const onTargetRender = spatialPerfEnabled
    ? ({ reason }: { reason: string }) => recordSpatialPerf(`densityGrid.targetRender.${reason}`, 0)
    : null;

/**
 * The layers every cell of a density grid draws, as one isolated group. The group draws them
 * into an offscreen target and composites that into each cell; cells under one view state, of
 * the same size to the device pixel, share the target, so the layers are drawn once a frame
 * however many cells are in view, rather than once per cell.
 *
 * Keep out of it what must reach its own canvas: the editable selection layer, which listens
 * on each canvas, and anything drawn per cell. Picks inside the group name the group as
 * `info.layer`; match them by member (`pickingInfoMatchesLayer`).
 */
export function sharedCellLayersGroup(layers: readonly Layer[]): IsolatedGroupLayer {
    return new IsolatedGroupLayer({
        id: SHARED_CELL_LAYERS_ID,
        layers: [...layers],
        onTargetRender,
        ...forwardedEvents,
    });
}
