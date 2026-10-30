import BaseChart from "@/charts/BaseChart";
import "@/charts/DensityScatterPlot";
import "@/charts/TableChart";
import "@/charts/WGLScatterPlot";
import "@/react/components/DeckScatterReactWrapper";
import "@/react/components/VivMDVReact";
import {
    convertLegacyChartConfig,
    convertLegacyChartConfigNow,
    getLegacyConversionTarget,
    isLegacyChartType,
} from "@/charts/legacyCharts";
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

        test("converts a classic table, keeping its sort in the new form", () => {
            convertLegacyChartsSetting.set("on");
            const config = { type: "table_chart", sort: { columnId: "x", sortAsc: false } };
            expect(convertLegacyChartConfig(config)).toBe("table_chart");
            expect(config).toEqual({ type: "table_chart_react", sort: { columnId: "x", ascending: false } });
        });

        test("converts a classic table that is not sorted", () => {
            convertLegacyChartsSetting.set("on");
            const config = { type: "table_chart" };
            expect(convertLegacyChartConfig(config)).toBe("table_chart");
            expect(config).toEqual({ type: "table_chart_react" });
        });

        test("converts a background filter, keeping the selected categories", () => {
            convertLegacyChartsSetting.set("on");
            const config = {
                type: "wgl_scatter_plot",
                background_filter: { column: "cell_type", categories: ["T", "B"], category: "T" },
            };
            expect(convertLegacyChartConfig(config)).toBe("wgl_scatter_plot");
            expect(config).toMatchObject({
                type: "DeckContourScatter",
                background_filter: { column: "cell_type", category: ["T", "B"] },
            });
        });

        test("converts a background filter saved with a single category", () => {
            convertLegacyChartsSetting.set("on");
            const config = { type: "wgl_scatter_plot", background_filter: { column: "cell_type", category: "T" } };
            convertLegacyChartConfig(config);
            expect(config).toMatchObject({
                type: "DeckContourScatter",
                background_filter: { column: "cell_type", category: ["T"] },
            });
        });

        test("drops a background filter that was not selecting anything", () => {
            convertLegacyChartsSetting.set("on");
            for (const background_filter of [
                { column: "__none__" },
                { column: null, categories: [] },
                { column: "cell_type", categories: [] },
            ]) {
                const config: { type: string; background_filter?: unknown } = { type: "wgl_scatter_plot", background_filter };
                expect(convertLegacyChartConfig(config)).toBe("wgl_scatter_plot");
                expect(config.background_filter).toBeUndefined();
            }
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

test("the classic table is legacy", () => {
    expect(BaseChart.types.table_chart.legacy).toBe(true);
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

describe("converting a chart on request", () => {
    test("names the replacement for a convertible chart, whatever the setting", () => {
        convertLegacyChartsSetting.set("off");
        expect(getLegacyConversionTarget({ type: "wgl_scatter_plot" })).toBe("DeckContourScatter");
        expect(getLegacyConversionTarget({ type: "table_chart" })).toBe("table_chart_react");
    });

    test("has no replacement for current charts or ones that cannot be converted", () => {
        expect(getLegacyConversionTarget({ type: "DeckContourScatter" })).toBeUndefined();
        expect(getLegacyConversionTarget({ type: "VivMdvRegionReact" })).toBeUndefined();
        expect(getLegacyConversionTarget({ type: "density_scatter_plot", viv: {} } as never)).toBeUndefined();
    });

    test("converts the config even when the load-time setting is off", () => {
        convertLegacyChartsSetting.set("off");
        const config = { type: "table_chart", sort: { columnId: "x", sortAsc: true } };
        expect(convertLegacyChartConfig({ ...config })).toBeUndefined();
        expect(convertLegacyChartConfigNow(config)).toBe("table_chart");
        expect(config).toEqual({ type: "table_chart_react", sort: { columnId: "x", ascending: true } });
    });
});

test("converting a classic scatter keeps its point size on screen", () => {
    const config: { type: string; radius: number; size: number[]; course_radius?: number } = {
        type: "wgl_scatter_plot",
        radius: 2,
        size: [493, 462],
    };
    convertLegacyChartConfigNow(config);
    expect(config).toMatchObject({ type: "DeckContourScatter", radius: 3.3, course_radius: 10 });
});

test("converting a deck scatter does not touch its point size", () => {
    const config = { type: "wgl_scatter_plot_dev", radius: 7, size: [493, 462] };
    convertLegacyChartConfigNow(config);
    expect(config).toEqual({ type: "DeckContourScatter", radius: 7, size: [493, 462] });
});
