import { MdvDataFilterExtension } from "@/webgl/MdvDataFilterExtension";
import { MdvScatterplotLayer } from "@/webgl/MdvScatterplotLayer";
import { ScatterDensityExension, ScatterSquareExtension } from "@/webgl/ScatterDeckExtension";
import { editWGSL } from "@/webgl/wgslLayerExtension";
import type { Layer, LayerExtension } from "@deck.gl/core";
import { ScatterplotLayer } from "@deck.gl/layers";
import type { ShaderModule, ShaderPlugin } from "@luma.gl/shadertools";
import { describe, expect, it } from "vitest";

type Shaders = { source: string; modules: ShaderModule[]; plugins?: ShaderPlugin[]; inject?: Record<string, string> };

/** Layer.getShaders reads `context`, which deck's LayerManager normally assigns. */
function withDevice<L extends Layer>(layer: L, type: "webgl" | "webgpu") {
    return Object.assign(layer, { context: { device: { type }, defaultShaderModules: [] } });
}

function shadersFor(extensions: LayerExtension[], type: "webgl" | "webgpu"): Shaders {
    return withDevice(new MdvScatterplotLayer({ id: "test", data: [], extensions }), type).getShaders();
}

const stockSource: string = withDevice(new ScatterplotLayer({ id: "stock", data: [] }), "webgpu").getShaders().source;

describe("editWGSL", () => {
    const source = "a;\nb;\nc;";
    it("inserts before, after, or replaces an anchor", () => {
        expect(editWGSL(source, { label: "t", anchor: "b;", position: "before", code: "x;" })).toBe("a;\nx;\nb;\nc;");
        expect(editWGSL(source, { label: "t", anchor: "b;", position: "after", code: "x;" })).toBe("a;\nb;\nx;\nc;");
        expect(editWGSL(source, { label: "t", anchor: "b;", position: "replace", code: "x;" })).toBe("a;\nx;\nc;");
    });
    it("throws when the anchor is missing or ambiguous", () => {
        expect(() => editWGSL(source, { label: "t", anchor: "z;", position: "after", code: "" })).toThrow(/anchor/);
        expect(() => editWGSL("b;b;", { label: "t", anchor: "b;", position: "after", code: "" })).toThrow(/anchor/);
    });
});

// These pin our anchors against the installed deck.gl: if an upgrade moves them, this
// fails here rather than as a render-time error on WebGPU.
describe("MdvScatterplotLayer on WebGPU", () => {
    it("applies every extension's WGSL edit to deck's scatterplot source", () => {
        const { source, plugins } = shadersFor(
            [new MdvDataFilterExtension({ filterSize: 1 }), new ScatterSquareExtension(), new ScatterDensityExension()],
            "webgpu",
        );
        expect(source).not.toBe(stockSource);
        expect(source).toContain("if (!dataFilter_isVisible(filterValues))");
        expect(source).toContain("max(abs(varyings.unitPosition.x), abs(varyings.unitPosition.y))");
        expect(source).not.toContain("let distToCenter = length(varyings.unitPosition)");
        expect(source).toContain("ScatterDensityExension.opacity");
        expect(plugins?.map((p) => p.vertexInputs)).toEqual([{ filterValues: "f32" }]);
    });

    it("gives the data filter a WGSL module and no GLSL injections", () => {
        const { modules, inject } = shadersFor([new MdvDataFilterExtension({ filterSize: 1 })], "webgpu");
        const dataFilter = modules.find((m) => m.name === "dataFilter");
        expect(dataFilter?.source).toContain("var<uniform> dataFilter: dataFilterUniforms");
        expect(dataFilter?.vs).toBeUndefined();
        expect(dataFilter?.uniformTypes).toMatchObject({ enabled: "i32", min: "vec4<f32>", max: "vec4<f32>" });
        expect(inject).toBeUndefined();
    });

    it("rejects DataFilterExtension options it has no WGSL for", () => {
        expect(() => shadersFor([new MdvDataFilterExtension({ filterSize: 1, categorySize: 1 })], "webgpu")).toThrow(
            /filterSize: 1/,
        );
    });

    it("leaves the shader alone without extensions", () => {
        expect(shadersFor([], "webgpu").source).toBe(stockSource);
    });
});

describe("MdvScatterplotLayer on WebGL", () => {
    it("matches ScatterplotLayer, with the extensions' GLSL injections", () => {
        const { source, plugins, inject, modules } = shadersFor(
            [new MdvDataFilterExtension({ filterSize: 1 }), new ScatterSquareExtension()],
            "webgl",
        );
        expect(source).toBe(stockSource);
        expect(plugins).toBeUndefined();
        expect(inject?.["fs:#main-end"]).toContain("ScatterSquareExtension");
        expect(modules.find((m) => m.name === "dataFilter")?.vs).toContain("dataFilterUniforms");
    });
});
