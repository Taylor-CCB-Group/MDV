import type { Layer } from "@deck.gl/core";
import { DataFilterExtension, type DataFilterExtensionProps } from "@deck.gl/extensions";
import type { ShaderModule, ShaderPlugin } from "@luma.gl/shadertools";
import { type WGSLLayerExtension, editWGSL, isWebGPULayer } from "./wgslLayerExtension";

/**
 * WGSL twin of the `dataFilter` uniform block from `@deck.gl/extensions`. The field
 * order must match upstream's `uniformTypes`, because luma lays the buffer out from
 * that table, and the variable must be named after the module (`dataFilter`).
 */
const dataFilterWGSL = /* wgsl */ `\
struct dataFilterUniforms {
  useSoftMargin: i32,
  enabled: i32,
  transformSize: i32,
  transformColor: i32,
  min: vec4<f32>,
  softMin: vec4<f32>,
  softMax: vec4<f32>,
  max: vec4<f32>,
  min64High: vec4<f32>,
  max64High: vec4<f32>,
  categoryBitMask: vec4<u32>,
};

@group(0) @binding(auto) var<uniform> dataFilter: dataFilterUniforms;

fn dataFilter_isVisible(value: f32) -> bool {
  return dataFilter.enabled == 0 || (value >= dataFilter.min.x && value <= dataFilter.max.x);
}
`;

/** Supplies `filterValues` (the attribute upstream's extension already manages) to WGSL. */
const filterValuesPlugin: ShaderPlugin = {
    name: "mdvDataFilterValues",
    vertexInputs: { filterValues: "f32" },
};

/**
 * `DataFilterExtension` that also works on WebGPU, for the subset of it MDV uses.
 *
 * On WebGL this is exactly upstream's extension. On WebGPU, where upstream only has
 * GLSL, it keeps upstream's JS side (the `filterValues` attribute and the uniform
 * values) and swaps the shader side for a WGSL module plus an edit to the host layer's
 * WGSL: objects whose value is outside `filterRange` collapse to a degenerate position,
 * as upstream does on WebGL.
 *
 * WebGPU supports only `filterSize: 1`, with no categories, `fp64` or `countItems`.
 * `filterSoftRange`, `filterTransformSize` and `filterTransformColor` are ignored there:
 * an object is either fully shown or hidden.
 *
 * Expected to be retired when deck.gl's DataFilterExtension supports WebGPU. The host
 * layer must be an `MdvScatterplotLayer`.
 */
export class MdvDataFilterExtension extends DataFilterExtension implements WGSLLayerExtension {
    static extensionName = "MdvDataFilterExtension";
    wgslPlugins = [filterValuesPlugin];

    getShaders(this: Layer<DataFilterExtensionProps>, extension: this) {
        const shaders = super.getShaders(extension);
        if (!isWebGPULayer(this)) return shaders;
        const { filterSize, categorySize, fp64, countItems } = extension.opts;
        if (filterSize !== 1 || categorySize || fp64 || countItems) {
            throw new Error(
                "MdvDataFilterExtension on WebGPU supports only { filterSize: 1 } " +
                    "(no categorySize, fp64 or countItems)",
            );
        }
        // Reuse upstream's module for its uniformTypes and getUniforms, so the uniform
        // values stay in step with upstream; replace only the shader source.
        const [upstream]: ShaderModule[] = shaders.modules;
        const module: ShaderModule = {
            name: upstream.name,
            uniformTypes: upstream.uniformTypes,
            getUniforms: upstream.getUniforms,
            source: dataFilterWGSL,
        };
        return { modules: [module] };
    }

    transformWGSL(source: string) {
        return editWGSL(source, {
            label: "MdvDataFilterExtension",
            anchor: "clip_filterPosition(&varyings.position, geometry.worldPosition.xy);",
            position: "after",
            code: `\
  // MdvDataFilterExtension
  if (!dataFilter_isVisible(filterValues)) {
    varyings.position = vec4<f32>(0.0);
  }`,
        });
    }
}
