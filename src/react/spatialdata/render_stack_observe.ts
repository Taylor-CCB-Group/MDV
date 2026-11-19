import type { LayerChannelConfig } from "@spatialdata/avivatorish";
import { type RenderStack, type RenderStackEntry, walkRenderStackEntries } from "@spatialdata/layers";

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

function stableRevisionValue(value: unknown): unknown {
    if (Array.isArray(value)) return value.map(stableRevisionValue);
    if (!value || typeof value !== "object") return value;

    const record = value as Record<string, unknown>;
    return Object.fromEntries(
        Object.keys(record)
            .sort()
            .map((key) => [key, stableRevisionValue(record[key])]),
    );
}

function spatialPropsState(props: Record<string, unknown>): Record<string, unknown> {
    const state: Record<string, unknown> = {};
    for (const key of Object.keys(props).sort()) {
        if (key === "channels" || key === "vivLayerProps") continue;
        state[key] = stableRevisionValue(props[key]);
    }
    return state;
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
    } else {
        void entry.props.opacity;
        void entry.props.blendMode;
        void entry.props.isolate;
        void entry.children.length;
    }
}

export function touchRenderStack(stack: RenderStack | undefined) {
    if (!stack) return;
    void stack.entries.length;
    for (const entry of walkRenderStackEntries(stack.entries)) {
        touchRenderStackEntry(entry);
    }
}

export function renderStackSpatialRevision(stack: RenderStack | undefined): string {
    if (!stack) return "";
    const parts: string[] = [];
    // Groups are written as brackets around their children: moving a layer out of a
    // group, just past its end, keeps the order but changes the opacity it inherits.
    const visit = (entries: readonly RenderStackEntry[]) => {
        for (const entry of entries) {
            if (entry.kind === "group") {
                // Spatial configs carry their groups' visibility and opacity.
                parts.push(
                    `${entry.id}:${entry.visible}:${String(entry.props.opacity ?? "")}:${String(entry.props.isolate ?? "")}:${String(entry.props.blendMode ?? "")}(`,
                );
                visit(entry.children);
                parts.push(")");
                continue;
            }
            if (entry.kind !== "spatial") continue;
            parts.push(`${entry.id}:${entry.visible}:${String(entry.props.opacity ?? "")}`);
            parts.push(JSON.stringify(spatialPropsState(entry.props)));
            const channels = entry.props.channels;
            if (channels && typeof channels === "object") {
                parts.push(channelStateKey(channels as LayerChannelConfig));
            }
            const vivLayerProps = entry.props.vivLayerProps;
            if (vivLayerProps && typeof vivLayerProps === "object") {
                parts.push(JSON.stringify(vivLayerProps));
            }
        }
    };
    visit(stack.entries);
    return parts.join("\0");
}
