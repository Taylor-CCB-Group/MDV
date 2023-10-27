
// --- copied straight from Avivator's code::: with notes for MDV ---

const DEFAUlT_CHANNEL_STATE = {
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

// --- following how VivViewMDV _parseChannels() works ---
export type MdvVivChannelConfig = {
    name: string,
    color?: `#${string}` | [r: number, g: number, b: number],
    visible?: boolean,
    contrastLimits?: [min: number, max: number],
    domains?: [min: number, max: number],
}
export type VivConfig = {
    channels: MdvVivChannelConfig[]
}

export type ROI = {
    min_x: number,
    min_y: number,
    max_x: number,
    max_y: number,
}

