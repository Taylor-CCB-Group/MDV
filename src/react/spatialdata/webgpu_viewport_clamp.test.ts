import { Framebuffer, type RenderPassParameters } from "@luma.gl/core";
import { describe, expect, it, vi } from "vitest";
import { applyWebGPUViewportClamp, attachmentSize, clampPassParameters } from "./webgpu_viewport_clamp";

/** A luma Framebuffer as its WebGPU subclass holds it: attachments are views of a texture. */
function makeFramebuffer(width: number, height: number, baseMipLevel = 0, reported = [width, height]) {
    return Object.assign(Object.create(Framebuffer.prototype), {
        width: reported[0],
        height: reported[1],
        colorAttachments: [{ props: { baseMipLevel }, texture: { width, height } }],
        depthStencilAttachment: null,
    });
}

describe("clampPassParameters", () => {
    it("clamps deck's y = -1 viewport on a backing store rounded down a pixel", () => {
        // CSS 1047x480 at dpr 2; devicePixelContentBox gave 959 rows, deck assumed 960
        expect(clampPassParameters({ viewport: [0, -1, 2094, 960, 0, 1] }, 2094, 959).viewport).toEqual([
            0, 0, 2094, 959, 0, 1,
        ]);
    });

    it("clamps a first frame drawn before the canvas leaves its default 300x150", () => {
        expect(clampPassParameters({ viewport: [0, -85.5, 250, 235.5] }, 300, 150).viewport).toEqual([
            0, 0, 250, 150, 0, 1,
        ]);
    });

    it("leaves a viewport inside the attachment as it is", () => {
        expect(clampPassParameters({ viewport: [0, 3, 500, 471, 0, 1] }, 500, 474).viewport).toEqual([
            0, 3, 500, 471, 0, 1,
        ]);
    });

    it("keeps the depth range", () => {
        expect(clampPassParameters({ viewport: [-10, 0, 120, 50, 0.25, 0.75] }, 100, 50).viewport).toEqual([
            0, 0, 100, 50, 0.25, 0.75,
        ]);
    });

    it("gives an empty rect for a viewport wholly outside", () => {
        expect(clampPassParameters({ scissorRect: [200, 10, 50, 20] }, 100, 50).scissorRect).toEqual([100, 10, 0, 20]);
    });

    it("clamps the scissor rect and passes other parameters through", () => {
        const parameters: RenderPassParameters = { scissorRect: [0, -1, 2094, 960], stencilReference: 1 };
        expect(clampPassParameters(parameters, 2094, 959)).toEqual({
            scissorRect: [0, 0, 2094, 959],
            stencilReference: 1,
        });
    });

    it("returns parameters without a viewport or scissor rect unchanged", () => {
        const parameters: RenderPassParameters = { blendConstant: [0, 0, 0, 1] };
        expect(clampPassParameters(parameters, 1, 1)).toBe(parameters);
    });
});

describe("attachmentSize", () => {
    it("uses the mip level's size, not the framebuffer's (luma's mipmap passes report 1x1)", () => {
        expect(attachmentSize(makeFramebuffer(1024, 512, 1, [1, 1]))).toEqual([512, 256]);
        expect(attachmentSize(makeFramebuffer(1024, 512, 10, [1, 1]))).toEqual([1, 1]);
    });
});

describe("applyWebGPUViewportClamp", () => {
    function makeDevice(type: string, framebuffer: unknown) {
        const setParameters = vi.fn();
        const pass = { framebuffer, setParameters };
        const device = Object.assign(Object.create(null), { type, beginRenderPass: vi.fn(() => pass) });
        return { device, setParameters };
    }

    it("clamps viewports set on a WebGPU device's render passes", () => {
        const { device, setParameters } = makeDevice("webgpu", makeFramebuffer(2094, 959));
        applyWebGPUViewportClamp(device);
        device.beginRenderPass({}).setParameters({ viewport: [0, -1, 2094, 960, 0, 1] });
        expect(setParameters).toHaveBeenCalledWith({ viewport: [0, 0, 2094, 959, 0, 1] });
    });

    it("leaves WebGL devices alone", () => {
        const { device } = makeDevice("webgl", makeFramebuffer(10, 10));
        const beginRenderPass = device.beginRenderPass;
        applyWebGPUViewportClamp(device);
        expect(device.beginRenderPass).toBe(beginRenderPass);
    });
});
