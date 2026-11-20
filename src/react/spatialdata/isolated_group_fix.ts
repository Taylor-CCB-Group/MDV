import type { Layer, LayersList, UpdateParameters } from "@deck.gl/core";
import { IsolatedGroupLayer } from "@spatialdata/layers";

/**
 * Runtime fix for `IsolatedGroupLayer` (@spatialdata/layers 0.17): a layer whose blend mode
 * is toggled inside an isolated group swaps state with its own wrapper.
 *
 * Inside an isolated group, `renderStackToDeckLayers` draws a layer with a blend mode as a
 * one-member `IsolatedGroupLayer` that takes the layer's own id. The parent group scopes its
 * members as `<group>/<member id>`, so the layer and its wrapper get the same sublayer id.
 * deck (9.4, `LayerManager._transferLayerState`) matches sublayers by id without checking
 * their class, so toggling the mode hands the layer's state to its wrapper, a group that
 * then never runs `initializeState`, or the wrapper's state to the layer. The parent compares
 * only member ids, so it sees no change of contents either. It happened to draw correctly in
 * testing; this keeps it from depending on that.
 *
 * The fix gives such a wrapper a sublayer id of its own (`<scoped id>-blend`: `-` keeps
 * `entryForLabel` attributing its labels to the entry), and treats a member that changes
 * class as a change of contents. Picks are unaffected: the wrapper, innermost, records the
 * member's own id. Remove once sd.js scopes the wrapper's id itself.
 */

export const BLEND_WRAPPER_SUFFIX = "-blend";

function flatten(layers: LayersList, into: Layer[] = []): Layer[] {
    for (const item of layers) {
        if (Array.isArray(item)) flatten(item, into);
        else if (item && typeof item === "object") into.push(item);
    }
    return into;
}

function layerName(layer: Layer): unknown {
    return Reflect.get(layer.constructor, "layerName");
}

/** Members as the parent sees them: id and class, so a layer and its wrapper differ. */
function memberSignature(layers: LayersList | undefined): string {
    return flatten(layers ?? [])
        .map((layer) => `${layer.id}:${String(layerName(layer))}`)
        .join("\0");
}

/** A one-member wrapper `renderStackToDeckLayers` made for a blended layer, as scoped by its parent. */
function isBlendWrapper(sublayer: Layer, scope: string): sublayer is IsolatedGroupLayer {
    if (!(sublayer instanceof IsolatedGroupLayer)) return false;
    const members = flatten(sublayer.props.layers);
    return members.length === 1 && sublayer.id === `${scope}${members[0].id}`;
}

/** The prefix `IsolatedGroupLayer` scopes its members' ids with (`/` and `%` escaped). */
function memberScope(groupId: string): string {
    return `${groupId.replaceAll("%", "%25").replaceAll("/", "%2F")}/`;
}

let applied = false;

export function applyIsolatedGroupFix(): void {
    if (applied) return;
    applied = true;
    const proto = IsolatedGroupLayer.prototype;
    const renderLayers = proto.renderLayers;
    const updateState = proto.updateState;
    proto.renderLayers = function (this: IsolatedGroupLayer): LayersList {
        const scope = memberScope(this.props.id);
        return flatten(renderLayers.call(this)).map((sublayer) =>
            isBlendWrapper(sublayer, scope)
                ? sublayer.clone({ id: `${sublayer.id}${BLEND_WRAPPER_SUFFIX}` })
                : sublayer,
        );
    };
    proto.updateState = function (this: IsolatedGroupLayer, params: UpdateParameters<IsolatedGroupLayer>): void {
        updateState.call(this, params);
        if (memberSignature(params.props.layers) !== memberSignature(params.oldProps.layers)) {
            this._markContentsChanged();
        }
    };
}
