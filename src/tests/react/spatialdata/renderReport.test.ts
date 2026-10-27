import type { DeviceInfo } from "@luma.gl/core";
import { describe, expect, test } from "vitest";
import { describeGpu, formatRenderReport } from "@/react/spatialdata/render_report";

const base = {
    renderer: "WebGPU",
    gpu: "webgpu: apple metal-3 (apple, integrated)",
    userAgent: "TestBrowser/1.0",
    url: "http://localhost/project/5",
};

describe("formatRenderReport", () => {
    test("keeps stacks as lines, and pairs a cause with how it was reported", () => {
        const report = formatRenderReport({
            ...base,
            failures: [
                {
                    name: "image: morphology_focus",
                    failure: {
                        reason: "Texture size exceeded",
                        detail: "Model(image).draw failed",
                        source: "webgpu",
                        count: 3,
                        dismissed: false,
                    },
                },
            ],
            canvasErrors: [
                { message: "boom", stack: "Error: boom\n    at a (x.ts:1)\n    at b (y.ts:2)", count: 2, dismissed: false },
            ],
        });
        expect(report).toContain("Renderer: WebGPU");
        expect(report).toContain("Layer turned off: image: morphology_focus (webgpu, ×3)\n    Texture size exceeded");
        expect(report).toContain("  Reported as:\n    Model(image).draw failed");
        expect(report).toContain("Rendering error (×2):\n    boom\n  Stack:\n    Error: boom\n        at a (x.ts:1)\n        at b (y.ts:2)");
        expect(report).not.toContain("\\n");
    });

    test("says so when there is nothing to report", () => {
        expect(formatRenderReport({ ...base, gpu: null, failures: [], canvasErrors: [] })).toContain(
            "GPU: not known yet",
        );
    });
});

describe("describeGpu", () => {
    test("leaves out unknown parts", () => {
        const info = {
            type: "webgpu",
            vendor: "apple",
            renderer: "metal-3",
            version: "",
            gpu: "apple",
            gpuType: "unknown",
            gpuBackend: "metal",
            shadingLanguage: "wgsl",
            shadingLanguageVersion: 100,
        } satisfies DeviceInfo;
        expect(describeGpu(info)).toBe("webgpu: apple metal-3 (apple, metal)");
    });
});
