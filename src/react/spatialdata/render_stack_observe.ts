import type { RenderStack, RenderStackEntry } from "@spatialdata/layers";
import type { LayerChannelConfig } from "@spatialdata/avivatorish";

function channelStateKey(channels: LayerChannelConfig): string {
    return JSON.stringify({
        channelIds: channels.channelIds,
        colors: channels.colors,
        contrastLimits: channels.contrastLimits,
        channelsVisible: channels.channelsVisible,
        selections: channels.selections,
    });
}

function touchChannelConfig(channels: unknown) {
    if (!channels || typeof channels !== "object") return;
    const config = channels as LayerChannelConfig;
    void config.channelIds;
    void config.colors;
    void config.contrastLimits;
    void config.channelsVisible;
    void config.selections;
}

function touchVivLayerProps(vivLayerProps: unknown) {
    if (!vivLayerProps || typeof vivLayerProps !== "object") return;
    const props = vivLayerProps as Record<string, unknown>;
    void props.brightness;
    void props.contrast;
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
        touchVivLayerProps(props.vivLayerProps);
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
            parts.push(channelStateKey(channels as LayerChannelConfig));
        }
        const vivLayerProps = entry.props.vivLayerProps;
        if (vivLayerProps && typeof vivLayerProps === "object") {
            parts.push(JSON.stringify(vivLayerProps));
        }
    }
    return parts.join("\0");
}
