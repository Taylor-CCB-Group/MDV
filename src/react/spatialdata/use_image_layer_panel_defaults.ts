import type { SpatialData } from "@spatialdata/core";
import { useMemo } from "react";

import {
    channelNamesFromImageChannels,
    loaderDefaultsFromImageChannels,
    type LoaderDefaults,
} from "@/react/spatialdata/image_layer_channel_bridge";

export function useImageLayerPanelDefaults(
    spatialData: SpatialData | null | undefined,
    elementKey: string | undefined,
): { channelNames: string[]; loaderDefaults: LoaderDefaults | undefined } {
    return useMemo(() => {
        if (!spatialData || !elementKey) {
            return { channelNames: [], loaderDefaults: undefined };
        }
        const image = spatialData.images?.[elementKey];
        if (!image) {
            return { channelNames: [], loaderDefaults: undefined };
        }
        const omeroChannels = image.channels;
        return {
            channelNames: channelNamesFromImageChannels(omeroChannels),
            loaderDefaults: loaderDefaultsFromImageChannels(omeroChannels),
        };
    }, [elementKey, spatialData]);
}
