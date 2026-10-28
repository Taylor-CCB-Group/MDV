import BaseChart from "@/charts/BaseChart";
import "@/charts/DensityScatterPlot";
import "@/charts/WGLScatterPlot";
import "@/react/components/DeckScatterReactWrapper";
import "@/react/components/VivMDVReact";
import { convertLegacyChartConfig, isLegacyChartType } from "@/charts/legacyCharts";
import { convertLegacyChartsSetting } from "@/lib/browserSettings";
import type DataStore from "@/datastore/DataStore";
import { afterEach, describe, expect, test } from "vitest";

afterEach(() => {
    convertLegacyChartsSetting.set("off");
});

function dataStoreWithRegions(regions: Record<string, { viv_image?: boolean; spatial?: boolean }>) {
    // only `regions.all_regions` is read by the chart type predicates
    return { regions: { all_regions: regions } } as unknown as DataStore;
}

describe("convertLegacyChartConfig", () => {
    test("leaves charts alone unless the setting is on", () => {
        const config = { type: "wgl_scatter_plot" };
        expect(convertLegacyChartConfig(config)).toBeUndefined();
        expect(config.type).toBe("wgl_scatter_plot");
    });

    describe("with the setting on", () => {
        test.each(["wgl_scatter_plot", "wgl_scatter_plot_dev", "density_scatter_plot"])(
            "%s becomes DeckContourScatter",
            (type) => {
                convertLegacyChartsSetting.set("on");
                const config = { type };
                expect(convertLegacyChartConfig(config)).toBe(type);
                expect(config.type).toBe("DeckContourScatter");
            },
        );

        test("does not convert a chart using a background filter the deck chart ignores", () => {
            convertLegacyChartsSetting.set("on");
            const config = { type: "wgl_scatter_plot", background_filter: { column: "cell_type", category: "T" } };
            expect(convertLegacyChartConfig(config)).toBeUndefined();
            expect(config.type).toBe("wgl_scatter_plot");
        });

        test("converts a chart whose background filter is the unset default", () => {
            convertLegacyChartsSetting.set("on");
            const config = { type: "wgl_scatter_plot", background_filter: { column: "__none__" } };
            expect(convertLegacyChartConfig(config)).toBe("wgl_scatter_plot");
        });

        test("leaves types with no equivalent alone", () => {
            convertLegacyChartsSetting.set("on");
            for (const type of ["VivMdvRegionReact", "wgl_3d_scatter_plot", "DeckContourScatter"]) {
                const config = { type };
                expect(convertLegacyChartConfig(config)).toBeUndefined();
                expect(config.type).toBe(type);
            }
        });
    });
});

describe("the Spatial Image Viewer", () => {
    const viewer = BaseChart.types.VivMdvRegionReact;

    test("is legacy when every region with a Viv image is also SpatialData", () => {
        const ds = dataStoreWithRegions({ a: { viv_image: true, spatial: true }, b: { spatial: true } });
        expect(isLegacyChartType(viewer, ds)).toBe(true);
    });

    test("is still offered for a Viv image region that is not SpatialData", () => {
        const ds = dataStoreWithRegions({ a: { viv_image: true, spatial: true }, b: { viv_image: true } });
        expect(isLegacyChartType(viewer, ds)).toBe(false);
    });
});
