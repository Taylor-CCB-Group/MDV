import { type PlatformInfo, type ShaderModule, WGSLShaderAssembler } from "@luma.gl/shadertools";
import { describe, expect, it, vi } from "vitest";
import { installWGSLAssemblyCache } from "./wgsl_assembly_cache";

const source = /* wgsl */ `
@vertex fn vertexMain(@location(0) position: vec2<f32>) -> @builtin(position) vec4<f32> {
  return vec4<f32>(tint(position), 0.0, 1.0);
}
@fragment fn fragmentMain() -> @location(0) vec4<f32> {
  return vec4<f32>(1.0);
}
`;

const tintModule: ShaderModule = {
    name: "tint",
    source: "fn tint(p: vec2<f32>) -> vec2<f32> { return p; }",
};

function platformInfo(): PlatformInfo {
    return {
        type: "webgpu",
        shaderLanguage: "wgsl",
        shaderLanguageVersion: 300,
        gpu: "apple",
        limits: {},
        features: new Set(),
    };
}

function cachedAssembler() {
    const assembler = new WGSLShaderAssembler();
    const assemble = vi.spyOn(assembler, "assembleWGSLShader");
    installWGSLAssemblyCache(assembler);
    return { assembler, assemble };
}

describe("installWGSLAssemblyCache", () => {
    it("assembles a repeated shader once and returns the same source", () => {
        const { assembler, assemble } = cachedAssembler();
        const device = platformInfo();
        const first = assembler.assembleWGSLShader({ platformInfo: device, source, modules: [tintModule] });
        const second = assembler.assembleWGSLShader({ platformInfo: device, source, modules: [tintModule] });
        expect(assemble).toHaveBeenCalledTimes(1);
        expect(second.source).toBe(first.source);
        expect(second.shaderLayout).toEqual(first.shaderLayout);
    });

    it("hands each caller its own layout objects", () => {
        const { assembler } = cachedAssembler();
        const device = platformInfo();
        const first = assembler.assembleWGSLShader({ platformInfo: device, source, modules: [tintModule] });
        const second = assembler.assembleWGSLShader({ platformInfo: device, source, modules: [tintModule] });
        expect(second.shaderLayout).not.toBe(first.shaderLayout);
        expect(second.bindingTable).not.toBe(first.bindingTable);
        expect(second.modules).not.toBe(first.modules);
    });

    it("assembles again when an input differs", () => {
        const { assembler, assemble } = cachedAssembler();
        const device = platformInfo();
        assembler.assembleWGSLShader({ platformInfo: device, source, modules: [tintModule] });
        assembler.assembleWGSLShader({ platformInfo: device, source, modules: [tintModule], defines: { FLAG: true } });
        assembler.assembleWGSLShader({ platformInfo: device, source, modules: [{ ...tintModule }] });
        // another device: same values, different `limits` / `features` objects
        assembler.assembleWGSLShader({ platformInfo: platformInfo(), source, modules: [tintModule] });
        expect(assemble).toHaveBeenCalledTimes(4);
    });
});
