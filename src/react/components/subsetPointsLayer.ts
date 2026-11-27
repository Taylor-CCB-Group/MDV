import type { Layer, PickingInfo } from "@deck.gl/core";
import SpatialLayer from "@/webgl/SpatialLayer";
import type { CategoryRows } from "./colorGridUtils";

/**
 * A points layer drawing only some of its rows: the layer with `rows` as its data, and only the
 * highlighted rows among them. Picks name a position in `rows`; the layer's handlers expect one in
 * all of its rows, so they are given that. `id`, when given, replaces the layer's (a grid cell's
 * layer takes the cell's view id as its prefix, to draw only there).
 */
export function subsetPointsLayer(
    points: Layer,
    { rows, positions }: CategoryRows,
    highlightedData: number[],
    id?: string,
): Layer {
    const toAllRows = (info: PickingInfo): PickingInfo =>
        info.index < 0 ? info : { ...info, index: positions[info.index] ?? -1 };
    const { onHover, onClick } = points.props;
    return points.clone({
        ...(id !== undefined && { id }),
        data: rows,
        // the highlight is drawn apart from `data`, so it has to be narrowed as well
        ...(points instanceof SpatialLayer && { highlightedData }),
        onHover: onHover
            ? (info: PickingInfo, event: Parameters<typeof onHover>[1]) => onHover(toAllRows(info), event)
            : null,
        onClick: onClick
            ? (info: PickingInfo, event: Parameters<typeof onClick>[1]) => onClick(toAllRows(info), event)
            : null,
    });
}
