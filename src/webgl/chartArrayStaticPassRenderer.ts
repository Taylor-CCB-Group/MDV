import { Deck, OrthographicView } from "@deck.gl/core";
import type { Layer, OrthographicViewState } from "@deck.gl/core";
import type { Device, Framebuffer } from "@luma.gl/core";
import { CHART_ARRAY_STATIC_VIEW_FLIP_Y } from "./chartArrayStaticBitmapLayers";

type RenderStaticPassArgs = {
    device: Device;
    target: Framebuffer;
    cssWidth: number;
    cssHeight: number;
    viewState: OrthographicViewState;
    layers: Layer[];
};

/**
 * Offscreen deck used to render chart-shared geometry into a framebuffer once.
 * Uses public deck props (`device`, `_framebuffer`) instead of internal methods.
 */
export class ChartArrayStaticPassRenderer {
    private deck: Deck<any> | null = null;
    private canvas: HTMLCanvasElement | null = null;
    private staticView = new OrthographicView({
        id: "chart-array-static-pass-view",
        controller: false,
        flipY: CHART_ARRAY_STATIC_VIEW_FLIP_Y,
    });

    render(args: RenderStaticPassArgs): void {
        if (typeof document === "undefined") return;
        if (!this.deck) {
            this.canvas = document.createElement("canvas");
            this.deck = new Deck<any>({
                id: "chart-array-static-pass-deck",
                canvas: this.canvas,
                width: args.cssWidth,
                height: args.cssHeight,
                useDevicePixels: true,
                controller: false,
                views: [this.staticView],
                viewState: {
                    ...args.viewState,
                } as any,
                layers: [],
                device: args.device,
                onError: null,
            });
        }

        this.deck?.setProps({
            width: args.cssWidth,
            height: args.cssHeight,
            useDevicePixels: true,
            views: [this.staticView],
            viewState: { ...args.viewState } as any,
            layers: args.layers,
            _framebuffer: args.target,
        } as any);
        this.deck?.redraw("chart-array-static-pass");
    }

    finalize() {
        this.deck?.finalize();
        this.deck = null;
        this.canvas = null;
    }
}

