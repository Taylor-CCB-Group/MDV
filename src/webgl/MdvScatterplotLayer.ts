import { ScatterplotLayer } from "@deck.gl/layers";
import type { ShaderPlugin } from "@luma.gl/shadertools";
import { isWGSLLayerExtension, isWebGPULayer } from "./wgslLayerExtension";

/**
 * deck's ScatterplotLayer, plus WGSL edits from our extensions when on WebGPU.
 *
 * On WebGL this behaves exactly like ScatterplotLayer: extensions inject GLSL through
 * deck's shader hooks. On WebGPU deck has no WGSL hooks, so each extension that
 * implements `WGSLLayerExtension` edits the WGSL source here instead.
 *
 * Use this (rather than ScatterplotLayer) for any scatterplot that carries MDV's
 * extensions, so that it can render on either backend.
 */
export class MdvScatterplotLayer<DataT = unknown, ExtraProps extends object = object> extends ScatterplotLayer<
    DataT,
    ExtraProps
> {
    static layerName = "MdvScatterplotLayer";

    getShaders() {
        const shaders = super.getShaders();
        if (!isWebGPULayer(this)) return shaders;
        let source: string = shaders.source;
        const plugins: ShaderPlugin[] = [...(shaders.plugins ?? [])];
        for (const extension of this.props.extensions) {
            if (!isWGSLLayerExtension(extension)) continue;
            source = extension.transformWGSL(source);
            if (extension.wgslPlugins) plugins.push(...extension.wgslPlugins);
        }
        return { ...shaders, source, plugins };
    }
}
