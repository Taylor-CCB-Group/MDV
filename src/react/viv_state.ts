import { useMemo } from "react";
import { useConfig } from "./hooks";
import { ColorPaletteExtension, loadOmeTiff } from "@hms-dbmi/viv";
import { useOmeTiff } from "./context";

export type OME_TIFF = Awaited<ReturnType<typeof loadOmeTiff>>;


// --- copied straight from Avivator's code::: with notes for MDV ---

export const DEFAUlT_CHANNEL_STATE = {
    channelsVisible: [],
    contrastLimits: [],
    colors: [],
    domains: [],
    selections: [],
    ids: [],
    // not for serialization... think about this.
    loader: [{ labels: [], shape: [] }],
    image: 0
};
const DEFAUlT_CHANNEL_VALUES = {
    channelsVisible: true,
    contrastLimits: [0, 65535],
    colors: [255, 255, 255],
    domains: [0, 65535],
    selections: { z: 0, c: 0, t: 0 },
    ids: ''
};

// --- following how VivViewerMDV _parseChannels() works ---
export type MdvVivChannelConfig = {
    name: string,
    color?: `#${string}` | [r: number, g: number, b: number],
    visible?: boolean,
    contrastLimits?: [min: number, max: number],
    domains?: [min: number, max: number],
}
// mobx for everything? (may make sense but I'm not a huge fan)
// we could have a config.channelsState, in which case it shouldn't need
// too much more architecture to make it work.
// Color
export type ChannelsState = {
    channelsVisible: boolean[],
    contrastLimits: [number, number][],
    colors: [number, number, number][],
    domains: [number, number][],
    selections: { z: number, c: number, t: number }[],
    ids: string[]
}

export type VivConfig = {
    channels: MdvVivChannelConfig[],
    image_properties:  ChannelsState,
    url: string
}

export type ROI = {
    min_x: number,
    min_y: number,
    max_x: number,
    max_y: number,
}

export function useChannelsState() {
    const config = useConfig<any>();
    //todo something useful... some values...
    return (config.image_properties || config.viv.image_properties) as ChannelsState;
}

export function useVivLayerConfig() {
    const imageProps = useChannelsState();
    const extensions = useMemo(() => [new ColorPaletteExtension()], []);

    const ome = useOmeTiff(); 
    if (!ome) return undefined;
    return {
        ...imageProps,
        loader: ome.data,
        extensions
    }
}
