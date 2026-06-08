import type { Device, Framebuffer, Texture } from "@luma.gl/core";

type ChartArrayStaticCacheKeyInput = {
    width: number;
    height: number;
    devicePixelRatio: number;
};

/** Buffer allocation key (size only — content is re-rendered in-place when this matches). */
export function getChartArrayStaticCacheKey({
    width,
    height,
    devicePixelRatio,
}: Pick<ChartArrayStaticCacheKeyInput, "width" | "height" | "devicePixelRatio">): string {
    const roundedWidth = Math.max(1, Math.round(width));
    const roundedHeight = Math.max(1, Math.round(height));
    const roundedDpr = Math.max(1, Math.round(devicePixelRatio * 1000) / 1000);
    return `${roundedWidth}x${roundedHeight}@${roundedDpr}`;
}

/** Content key — pan/zoom and layer set; triggers FBO repaint without reallocating the texture. */
export function getChartArrayStaticContentKey(
    viewState: { target?: number[]; zoom?: number },
    staticLayerIds: readonly string[],
): string {
    const target = viewState.target ?? [0, 0, 0];
    const roundedTarget = target
        .map((value) => Math.round(Number(value) * 1000) / 1000)
        .join(",");
    const zoomValue = viewState.zoom;
    const roundedZoom = Array.isArray(zoomValue)
        ? zoomValue.map((value) => Math.round(Number(value) * 1000) / 1000).join(",")
        : String(
              Math.round((Number.isFinite(zoomValue) ? Number(zoomValue) : 0) * 1000) / 1000,
          );
    const ids = [...staticLayerIds].sort().join("\u0001");
    return `${roundedTarget}|z=${roundedZoom}|${ids}`;
}

export type ChartArrayStaticPass = {
    cacheKey: string;
    /** Physical pixel width of the color attachment. */
    width: number;
    /** Physical pixel height of the color attachment. */
    height: number;
    /** CSS pixel width used for orthographic projection (matches grid cells). */
    cssWidth: number;
    /** CSS pixel height used for orthographic projection (matches grid cells). */
    cssHeight: number;
    framebuffer: Framebuffer;
    colorTexture: Texture;
    reused: boolean;
};

type ChartArrayRenderCacheEntry = {
    cacheKey: string;
    width: number;
    height: number;
    cssWidth: number;
    cssHeight: number;
    framebuffer: Framebuffer;
    colorTexture: Texture;
};

export class ChartArrayRenderCache {
    private staticPass: ChartArrayRenderCacheEntry | null = null;

    getOrCreateStaticPass(args: {
        device: Device;
        cacheKey: string;
        /** Physical pixel width of the color attachment. */
        width: number;
        /** Physical pixel height of the color attachment. */
        height: number;
        /** CSS pixel width for orthographic projection (matches grid cells). */
        cssWidth: number;
        /** CSS pixel height for orthographic projection (matches grid cells). */
        cssHeight: number;
    }): ChartArrayStaticPass {
        const width = Math.max(1, Math.round(args.width));
        const height = Math.max(1, Math.round(args.height));
        const matchesExisting =
            this.staticPass &&
            this.staticPass.cacheKey === args.cacheKey &&
            this.staticPass.width === width &&
            this.staticPass.height === height;

        if (matchesExisting && this.staticPass) {
            return {
                ...this.staticPass,
                reused: true,
            };
        }

        this.dispose();

        const colorTexture = args.device.createTexture({
            format: "rgba8unorm",
            width,
            height,
        });
        const framebuffer = args.device.createFramebuffer({
            width,
            height,
            colorAttachments: [colorTexture],
        });
        const cssWidth = Math.max(1, Math.round(args.cssWidth ?? width));
        const cssHeight = Math.max(1, Math.round(args.cssHeight ?? height));
        this.staticPass = {
            cacheKey: args.cacheKey,
            width,
            height,
            cssWidth,
            cssHeight,
            framebuffer,
            colorTexture,
        };
        return {
            ...this.staticPass,
            reused: false,
        };
    }

    invalidateStaticPass() {
        this.dispose();
    }

    dispose() {
        if (!this.staticPass) return;
        this.staticPass.framebuffer.destroy();
        this.staticPass.colorTexture.destroy();
        this.staticPass = null;
    }
}
