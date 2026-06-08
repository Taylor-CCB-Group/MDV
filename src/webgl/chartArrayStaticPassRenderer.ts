import { Deck, OrthographicView } from "@deck.gl/core";
import type { Layer, OrthographicViewState } from "@deck.gl/core";
import type { Device, Framebuffer } from "@luma.gl/core";
import { CHART_ARRAY_STATIC_VIEW_FLIP_Y } from "./chartArrayStaticBitmapLayers";

const STATIC_PASS_VIEW_ID = "chart-array-static-pass-view";

type RenderStaticPassArgs = {
    device: Device;
    target: Framebuffer;
    cssWidth: number;
    cssHeight: number;
    viewState: OrthographicViewState;
    layers: Layer[];
};

function createStaticPassView(cssWidth: number, cssHeight: number): OrthographicView {
    return new OrthographicView({
        id: STATIC_PASS_VIEW_ID,
        x: 0,
        y: 0,
        width: cssWidth,
        height: cssHeight,
        controller: false,
        flipY: CHART_ARRAY_STATIC_VIEW_FLIP_Y,
    });
}

/**
 * Offscreen deck used to render chart-shared geometry into a framebuffer once.
 * Uses public deck props (`device`, `_framebuffer`) instead of internal methods.
 */
export class ChartArrayStaticPassRenderer {
    private deck: Deck<any> | null = null;
    private canvas: HTMLCanvasElement | null = null;

    render(args: RenderStaticPassArgs): void {
        if (typeof document === "undefined") return;
        const staticView = createStaticPassView(args.cssWidth, args.cssHeight);
        const viewState = {
            ...args.viewState,
            id: STATIC_PASS_VIEW_ID,
        };
        if (!this.deck) {
            this.canvas = document.createElement("canvas");
            this.deck = new Deck<any>({
                id: "chart-array-static-pass-deck",
                canvas: this.canvas,
                width: args.cssWidth,
                height: args.cssHeight,
                useDevicePixels: true,
                controller: false,
                views: [staticView],
                viewState,
                layers: [],
                device: args.device,
                onError: null,
            });
        }

        this.deck?.setProps({
            width: args.cssWidth,
            height: args.cssHeight,
            useDevicePixels: true,
            views: [staticView],
            viewState,
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
