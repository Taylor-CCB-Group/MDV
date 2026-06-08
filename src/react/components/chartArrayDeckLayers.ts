import type { Layer } from "@deck.gl/core";
import {
    getDeckLayerViewportScope,
    tagDeckLayerViewportScope,
} from "./deckLayerViewportScope";
import { cloneLayerForChartArrayGrid } from "./densityGridUtils";
import {
    composeLayerGroups,
    type LayerGroup,
} from "@/webgl/layerGroups";

type CloneableDeckLayer = Layer & {
    clone: (props: Record<string, unknown>) => Layer;
};

export type BuildChartArrayDeckLayersOptions = {
    visibleCellIndices: readonly number[];
    viewIds: readonly string[];
    greyScatterplotLayer: CloneableDeckLayer | null;
    scatterplotLayer: CloneableDeckLayer | null;
    chartSharedLayers: readonly (Layer | null)[];
    buildPerViewportLayer: (
        index: number,
        viewId: string,
    ) => Layer | null;
};

function asChartSharedLayer(layer: CloneableDeckLayer): Layer {
    if (getDeckLayerViewportScope(layer)) {
        return layer;
    }
    return tagDeckLayerViewportScope(layer, "chart-shared");
}

function cloneSharedLayerForGrid(
    layer: CloneableDeckLayer,
    props: Record<string, unknown> = {},
): Layer {
    return asChartSharedLayer(cloneLayerForChartArrayGrid(layer, props) as CloneableDeckLayer);
}

function adaptChartSharedLayerForGrid(layer: Layer): Layer {
    const cloneable = layer as CloneableDeckLayer;
    if (!cloneable.clone) return layer;
    return cloneSharedLayerForGrid(cloneable);
}

export function buildChartArrayDeckLayers({
    visibleCellIndices,
    viewIds,
    greyScatterplotLayer,
    scatterplotLayer,
    chartSharedLayers,
    buildPerViewportLayer,
}: BuildChartArrayDeckLayersOptions): Layer[] {
    const groups: LayerGroup[] = [
        {
            id: "shared-geometry",
            scope: "chart-shared",
            buildLayers: () => {
                const sharedGeometryLayers: Layer[] = [];
                if (greyScatterplotLayer) {
                    sharedGeometryLayers.push(cloneSharedLayerForGrid(greyScatterplotLayer));
                }
                if (scatterplotLayer) {
                    sharedGeometryLayers.push(
                        cloneSharedLayerForGrid(scatterplotLayer, { contourLayers: [] }),
                    );
                }
                return sharedGeometryLayers;
            },
        },
        {
            id: "density-per-viewport",
            scope: "per-viewport",
            buildLayers: (context) =>
                context.visibleCellIndices.flatMap((index) => {
                    const viewId = context.viewIds[index];
                    if (!viewId) return [];
                    const layer = buildPerViewportLayer(index, viewId);
                    return layer ? [layer] : [];
                }),
        },
        {
            id: "chart-shared-overlays",
            scope: "chart-shared",
            buildLayers: () =>
                chartSharedLayers
                    .filter((layer): layer is Layer => layer !== null)
                    .map(adaptChartSharedLayerForGrid),
        },
    ];
    return composeLayerGroups(groups, { visibleCellIndices, viewIds });
}
