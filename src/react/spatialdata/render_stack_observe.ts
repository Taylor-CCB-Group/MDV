import type { RenderStack, RenderStackEntry } from "@spatialdata/layers";

import {
    type ChannelConfig,
    channelStateKey,
} from "./image_layer_channel_bridge";

function touchChannelConfig(channels: unknown) {
    if (!channels || typeof channels !== "object") return;
    const config = channels as ChannelConfig;
    void config.channelIds;
    void config.colors;
    void config.contrastLimits;
    void config.channelsVisible;
    void config.selections;
}

/** Establish MobX subscriptions for a single stack entry (call during observer render). */
export function touchRenderStackEntry(entry: RenderStackEntry | undefined) {
    if (!entry) return;
    void entry.id;
    void entry.visible;
    void entry.kind;
    if (entry.kind === "spatial") {
        void entry.source.elementType;
        void entry.source.elementKey;
        const { props } = entry;
        void props.opacity;
        touchChannelConfig(props.channels);
        for (const key of Object.keys(props)) {
            void props[key];
        }
    } else if (entry.kind === "host") {
        void entry.source.hostLayerId;
    }
}

export function touchRenderStack(stack: RenderStack | undefined) {
    if (!stack) return;
    void stack.entries.length;
    for (const entry of stack.entries) {
        touchRenderStackEntry(entry);
    }
}

export function renderStackSpatialRevision(stack: RenderStack | undefined): string {
    if (!stack) return "";
    const parts: string[] = [];
    for (const entry of stack.entries) {
        if (entry.kind !== "spatial") continue;
        parts.push(`${entry.id}:${entry.visible}:${String(entry.props.opacity ?? "")}`);
        const channels = entry.props.channels;
        if (channels && typeof channels === "object") {
            parts.push(channelStateKey(channels as ChannelConfig));
        }
    }
    return parts.join("\0");
}
