import type { Layer, LayersList } from "@deck.gl/core";

export type LayerGroupScope = "chart-shared" | "per-viewport";

export type LayerGroupBlendMode =
    | "normal"
    | "additive"
    | "multiply"
    | "screen"
    | "custom";

export type LayerGroupRenderTarget = "screen" | "framebuffer";

export type LayerGroupBuildContext = {
    visibleCellIndices: readonly number[];
    viewIds: readonly string[];
};

export type LayerGroup = {
    id: string;
    scope: LayerGroupScope;
    blendMode?: LayerGroupBlendMode;
    renderTarget?: LayerGroupRenderTarget;
    buildLayers: (context: LayerGroupBuildContext) => LayersList;
};

export function composeLayerGroups(
    groups: readonly LayerGroup[],
    context: LayerGroupBuildContext,
): Layer[] {
    const layers: Layer[] = [];
    const appendLayers = (input: LayersList) => {
        for (const layer of input) {
            if (!layer) continue;
            if (Array.isArray(layer)) {
                appendLayers(layer);
            } else {
                layers.push(layer);
            }
        }
    };
    for (const group of groups) {
        appendLayers(group.buildLayers(context));
    }
    return layers;
}
