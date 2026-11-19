import { type RenderStack, walkRenderStackEntries } from "@spatialdata/layers";
import type { ImageLoaderData, LayerLoadState } from "@spatialdata/vis";

export type ImageLayerRegistry = {
    getImageLoadedDataByElementKey: (elementKey: string) => ImageLoaderData | undefined;
    getLayerLoadStateByElementKey?: (elementKey: string) => LayerLoadState | undefined;
};

export function getLayerLoadStateByElementKey(
    elementKey: string,
    stack: RenderStack,
    getLayerLoadState: (layerId?: string) => LayerLoadState | undefined,
): LayerLoadState | undefined {
    // The same element can be in several entries, some never loaded (e.g. in a hidden
    // group): the first one with a load state answers.
    for (const entry of walkRenderStackEntries(stack.entries)) {
        if (
            entry.kind === "spatial" &&
            entry.source.elementType === "image" &&
            entry.source.elementKey === elementKey
        ) {
            const state = getLayerLoadState(entry.id);
            if (state) return state;
        }
    }
    return undefined;
}

export function createImageLayerRegistry(
    stack: RenderStack,
    getImageLoadedDataByElementKey: ImageLayerRegistry["getImageLoadedDataByElementKey"],
    getLayerLoadState: (layerId?: string) => LayerLoadState | undefined,
): ImageLayerRegistry {
    return {
        getImageLoadedDataByElementKey,
        getLayerLoadStateByElementKey: (elementKey) =>
            getLayerLoadStateByElementKey(elementKey, stack, getLayerLoadState),
    };
}
