import { DetailView, ImageLayer, MultiscaleImageLayer } from "@hms-dbmi/viv";
import type { Layer } from "@deck.gl/core";

/** Synthetic DetailView id for the offscreen chart-array static raster pass. */
export const CHART_ARRAY_STATIC_PASS_VIEW_ID = "chart-array-static-pass-view";

export function isVivImageLayer(layer: Layer): boolean {
    return layer instanceof MultiscaleImageLayer || layer instanceof ImageLayer;
}

/** Layers rasterized once into the shared static FBO (viv image + scatter). */
export function isChartArrayStaticRasterLayer(layer: Layer): boolean {
    if (isVivImageLayer(layer)) {
        return true;
    }
    return layer.id.startsWith("scatter_") || layer.id.startsWith("scatter-grey_");
}

export function buildVivImageLayersForStaticPass(args: {
    width: number;
    height: number;
    layerConfig: Record<string, unknown>;
}): Layer[] {
    if (args.width <= 0 || args.height <= 0) {
        return [];
    }
    const view = new DetailView({
        id: CHART_ARRAY_STATIC_PASS_VIEW_ID,
        x: 0,
        y: 0,
        width: args.width,
        height: args.height,
    });
    return view.getLayers({ props: args.layerConfig }) as unknown as Layer[];
}
