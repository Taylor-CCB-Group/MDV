import type { SpatialData } from "@spatialdata/core";
import { useMemo } from "react";

import { useChart } from "@/react/context";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "@/react/components/SpatialDataMDVReact";
import {
    channelNamesFromImageChannels,
    loaderDefaultsFromImageChannels,
    type LoaderDefaults,
} from "@/react/spatialdata/image_layer_channel_bridge";

export function useImageLayerPanelDefaults(
    spatialData: SpatialData | null | undefined,
    elementKey: string | undefined,
) {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const registry = chart.imageLayerRegistry;

    return useMemo(() => {
        if (!elementKey) {
            return { channelNames: [] as string[], loaderDefaults: undefined as LoaderDefaults | undefined };
        }

        const registryData = registry?.getImageLoadedDataByElementKey(elementKey);
        if (registryData) {
            return {
                channelNames: registryData.channelNames ?? [],
                loaderDefaults: {
                    colors: registryData.colors,
                    contrastLimits: registryData.contrastLimits,
                    channelsVisible: registryData.channelsVisible,
                    selections: registryData.selections,
                } satisfies LoaderDefaults,
            };
        }

        const imageElement = spatialData?.images?.[elementKey] as
            | { channels?: { label?: string; color?: string; active?: boolean; window?: { start: number; end: number } }[] }
            | undefined;
        const omeroChannels = imageElement?.channels;
        if (omeroChannels?.length) {
            return {
                channelNames: channelNamesFromImageChannels(omeroChannels),
                loaderDefaults: loaderDefaultsFromImageChannels(omeroChannels),
            };
        }

        return { channelNames: [] as string[], loaderDefaults: undefined as LoaderDefaults | undefined };
    }, [elementKey, registry, spatialData]);
}
