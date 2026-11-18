import type { Layer } from "@deck.gl/core";
import type { ShaderPlugin } from "@luma.gl/shadertools";

/**
 * MDV-local stand-in for WGSL shader hooks, which deck.gl 9.4 does not have
 * (`SHADER_HOOKS_WGSL = [] // Not yet supported` in `@deck.gl/core` shaderlib).
 *
 * On WebGL our extensions inject GLSL through deck's hooks as usual. On WebGPU those
 * injections would be spliced verbatim into WGSL and fail to compile, so instead an
 * extension implementing this interface returns only WGSL shader modules from
 * `getShaders()`, and `MdvScatterplotLayer` asks it to edit the layer's WGSL source.
 *
 * Expected to be retired once deck ships WGSL hooks (and WGSL for `DataFilterExtension`);
 * the same pattern is used by `@spatialdata/layers` (`PointsScatterplotLayer`).
 */
export interface WGSLLayerExtension {
    /** Edit the host layer's WGSL source. Only called when the layer is on a WebGPU device. */
    transformWGSL(source: string): string;
    /** luma.gl shader plugins (e.g. extra vertex inputs) the edited source depends on. */
    wgslPlugins?: ShaderPlugin[];
}

export function isWGSLLayerExtension(extension: unknown): extension is WGSLLayerExtension {
    return (
        typeof extension === "object" &&
        extension !== null &&
        "transformWGSL" in extension &&
        typeof extension.transformWGSL === "function"
    );
}

export function isWebGPULayer(layer: Layer): boolean {
    return layer.context?.device?.type === "webgpu";
}

type WGSLEdit = {
    /** Used in the error message when the anchor is missing. */
    label: string;
    /** Exact text that must occur once in the source. */
    anchor: string;
    position: "before" | "after" | "replace";
    code: string;
};

/**
 * Apply an edit at an exact anchor in a deck.gl WGSL source.
 *
 * Throws rather than silently skipping when the anchor is absent or ambiguous: a deck
 * upgrade that moves the anchor should fail loudly (and `wgslLayerExtensions.test.ts`
 * pins the anchors against the installed deck), not render unfiltered or wrong-shaped
 * points.
 */
export function editWGSL(source: string, { label, anchor, position, code }: WGSLEdit): string {
    const index = source.indexOf(anchor);
    if (index < 0 || source.indexOf(anchor, index + anchor.length) >= 0) {
        throw new Error(
            `${label}: expected exactly one occurrence of WGSL anchor ${JSON.stringify(anchor)} ` +
                "(has the deck.gl shader changed?)",
        );
    }
    const end = index + anchor.length;
    switch (position) {
        case "before":
            return `${source.slice(0, index)}${code}\n${source.slice(index)}`;
        case "after":
            return `${source.slice(0, end)}\n${code}${source.slice(end)}`;
        case "replace":
            return `${source.slice(0, index)}${code}${source.slice(end)}`;
    }
}
