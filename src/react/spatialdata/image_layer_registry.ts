import type { RenderStack } from "@spatialdata/layers";
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
    for (const entry of stack.entries) {
        if (
            entry.kind === "spatial" &&
            entry.source.elementType === "image" &&
            entry.source.elementKey === elementKey
        ) {
            return getLayerLoadState(entry.id);
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
