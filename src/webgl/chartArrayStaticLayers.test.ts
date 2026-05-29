import { ImageLayer } from "@hms-dbmi/viv";
import { describe, expect, test } from "vitest";
import {
    CHART_ARRAY_STATIC_PASS_VIEW_ID,
    isChartArrayStaticRasterLayer,
    isVivImageLayer,
} from "./chartArrayStaticLayers";

describe("chartArrayStaticLayers", () => {
    test("identifies scatter layers as static raster layers", () => {
        expect(isChartArrayStaticRasterLayer({ id: "scatter_chart-array-grid" } as never)).toBe(true);
        expect(isChartArrayStaticRasterLayer({ id: "scatter-grey_chart-array-grid" } as never)).toBe(
            true,
        );
        expect(isChartArrayStaticRasterLayer({ id: "gate_-#detail#" } as never)).toBe(false);
    });

    test("identifies viv image layers by deck layer class", () => {
        const layer = Object.create(ImageLayer.prototype) as InstanceType<typeof ImageLayer>;
        layer.id = "ZarrImageLayer-#test#";
        expect(isVivImageLayer(layer)).toBe(true);
        expect(isChartArrayStaticRasterLayer(layer)).toBe(true);
    });

    test("uses a dedicated static-pass view id", () => {
        expect(CHART_ARRAY_STATIC_PASS_VIEW_ID).toBe("chart-array-static-pass-view");
    });
});
