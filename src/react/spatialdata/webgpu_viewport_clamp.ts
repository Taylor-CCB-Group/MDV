import { type Device, Framebuffer, type RenderPassParameters } from "@luma.gl/core";

/**
 * Keep deck's render-pass viewport and scissor rect inside the pass's attachment on WebGPU.
 *
 * deck 9.4 places each view with `getGLViewport` (`@deck.gl/core` `passes/layers-pass.js`):
 * `y = drawingBufferHeight - (viewport.y + viewport.height) * pixelRatio`, with the viewport
 * from the canvas's CSS size. Whenever the backing store disagrees with CSS size × dpr, the
 * rect falls outside the attachment:
 * - on a chart's first frames, before luma has resized the canvas from its default
 *   300×150 (`setViewport(0, -85, 250, 236)` on a 300×150 target, seen in MDV);
 * - when the CSS height is fractional and `devicePixelContentBox` rounds the backing store
 *   down by a pixel (`setViewport(0, -1, 2094, 960)` on 2094×959).
 *
 * Chrome accepts that. Safari 26/27 rejects it, which invalidates the whole render pass, so
 * nothing in it draws, and reports it as "GPUCommandEncoder.finish: encoder state is
 * 'Locked'" (a WebKit bug hides the real validation error). MDV then took the canvas as
 * failed and fell back to WebGL.
 *
 * Clamping costs at most a pixel of vertical scale while the sizes disagree, and on the
 * first frames draws into a canvas that is about to be resized anyway. Remove once deck
 * clamps `getGLViewport` to the drawing buffer itself.
 *
 * Wraps `beginRenderPass` on this device instance only; a no-op on WebGL.
 */
export function applyWebGPUViewportClamp(device: Device): void {
    if (device.type !== "webgpu") return;
    const beginRenderPass = device.beginRenderPass.bind(device);
    device.beginRenderPass = (props) => {
        const pass = beginRenderPass(props);
        // luma's WebGPURenderPass resolves the attachment (the canvas's current framebuffer
        // when none is given) in its constructor; it is not part of the public RenderPass type.
        const framebuffer: unknown = Reflect.get(pass, "framebuffer");
        const size = framebuffer instanceof Framebuffer ? attachmentSize(framebuffer) : null;
        if (!size) return pass;
        const setParameters = pass.setParameters.bind(pass);
        pass.setParameters = (parameters) => setParameters(clampPassParameters(parameters, size[0], size[1]));
        return pass;
    };
}

/**
 * The size of the mip level a framebuffer renders to. Not `framebuffer.width/height`: a
 * framebuffer built from views without a size (luma's WebGPU mipmap generation) reports 1×1.
 */
export function attachmentSize(framebuffer: Framebuffer): [width: number, height: number] | null {
    const view = framebuffer.colorAttachments[0] ?? framebuffer.depthStencilAttachment;
    if (!view) return null;
    const level = view.props.baseMipLevel ?? 0;
    return [Math.max(1, view.texture.width >> level), Math.max(1, view.texture.height >> level)];
}

/** [x, y, width, height] intersected with [0, 0, width, height] of the attachment. */
function clampRect(x: number, y: number, w: number, h: number, width: number, height: number) {
    const x0 = Math.min(Math.max(x, 0), width);
    const y0 = Math.min(Math.max(y, 0), height);
    const x1 = Math.min(Math.max(x + w, x0), width);
    const y1 = Math.min(Math.max(y + h, y0), height);
    return [x0, y0, x1 - x0, y1 - y0] as const;
}

/** `parameters` with `viewport` and `scissorRect` clamped to a `width`×`height` attachment. */
export function clampPassParameters(
    parameters: RenderPassParameters,
    width: number,
    height: number,
): RenderPassParameters {
    const { viewport, scissorRect } = parameters;
    if (!viewport && !scissorRect) return parameters;
    const clamped = { ...parameters };
    if (viewport) {
        const [x, y, w, h] = clampRect(viewport[0], viewport[1], viewport[2], viewport[3], width, height);
        clamped.viewport = [x, y, w, h, viewport[4] ?? 0, viewport[5] ?? 1];
    }
    if (scissorRect) {
        const [x, y, w, h] = clampRect(scissorRect[0], scissorRect[1], scissorRect[2], scissorRect[3], width, height);
        clamped.scissorRect = [x, y, w, h];
    }
    return clamped;
}
