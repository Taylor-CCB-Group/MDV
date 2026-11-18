import BaseChart from "@/charts/BaseChart";
import "@/charts/DensityScatterPlot";
import "@/charts/TableChart";
import "@/charts/WGLScatterPlot";
import "@/react/components/DeckScatterReactWrapper";
import "@/react/components/VivMDVReact";
import { afterEach, describe, expect, test } from "vitest";
import {
    convertLegacyChartConfig,
    convertLegacyChartConfigNow,
    getLegacyConversionTarget,
    isLegacyChartType,
} from "@/charts/legacyCharts";
import type DataStore from "@/datastore/DataStore";
import { convertLegacyChartsSetting } from "@/lib/browserSettings";

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
                const config: { type: string; background_filter?: unknown } = {
                    type: "wgl_scatter_plot",
                    background_filter,
                };
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

describe("the Viv scatter plot that is not a React chart", () => {
    const withRegions = (regions: Record<string, { viv_image?: unknown; spatial?: unknown }>) =>
        ({ regions: { all_regions: regions } }) as unknown as DataStore;

    test("converts for a region with a Viv image that is not SpatialData", () => {
        const dataStore = withRegions({ r1: { viv_image: {} } });
        expect(getLegacyConversionTarget({ type: "viv_scatter_plot", region: "r1" } as never, dataStore)).toBe(
            "VivMdvRegionReact",
        );
        const config = { type: "viv_scatter_plot", region: "r1", viv: { channels: [] } };
        expect(convertLegacyChartConfigNow(config, dataStore)).toBe("viv_scatter_plot");
        expect(config).toEqual({ type: "VivMdvRegionReact", region: "r1", viv: { channels: [] } });
    });

    test("converts a SpatialData region to the SpatialData.js viewer", () => {
        const dataStore = withRegions({ r1: { viv_image: {}, spatial: {} } });
        expect(getLegacyConversionTarget({ type: "viv_scatter_plot", region: "r1" } as never, dataStore)).toBe(
            "SpatialDataMdvRegionReact",
        );
        const config = { type: "viv_scatter_plot", region: "r1", radius: 0.5, viv: { channels: [] } };
        convertLegacyChartConfigNow(config, dataStore);
        expect(config).toEqual({ type: "SpatialDataMdvRegionReact", region: "r1", radius: 0.5 });
    });

    test("leaves a chart whose region has no Viv image, or that cannot be checked", () => {
        expect(
            getLegacyConversionTarget({ type: "viv_scatter_plot", region: "r1" } as never, withRegions({ r1: {} })),
        ).toBeUndefined();
        expect(getLegacyConversionTarget({ type: "viv_scatter_plot", region: "r1" } as never)).toBeUndefined();
    });
});

describe("the Spatial Image Viewer on a SpatialData region", () => {
    const withRegions = (regions: Record<string, { viv_image?: unknown; spatial?: unknown }>) =>
        ({ regions: { all_regions: regions } }) as unknown as DataStore;
    const spatialRegion = withRegions({
        r1: { viv_image: {}, spatial: { file: "a.zarr", coordinate_system: "global" } },
    });

    test("becomes the SpatialData.js viewer, keeping the scatter, region and filter settings", () => {
        const background_filter = { column: "sample", category: "r1" };
        const category_filters = [{ column: "cell_type", category: ["T"] }];
        const config = {
            type: "VivMdvRegionReact",
            id: "c1",
            title: "r1-image",
            region: "r1",
            param: ["x", "y"],
            color_by: "cell_type",
            radius: 4,
            opacity: 0.6,
            background_filter,
            category_filters,
            contourParameter: "cell_type",
            roi: { min_x: 0, max_x: 10, min_y: 0, max_y: 10 },
            json: "r1.json",
            showJson: true,
            viv: {
                channels: [{ name: "DAPI" }],
                viewerStore: { viewState: { target: [5, 5, 0], zoom: -2 } },
                channelsStore: { colors: [[255, 0, 0]], contrastLimits: [[0, 100]] },
            },
        };
        expect(convertLegacyChartConfigNow(config, spatialRegion)).toBe("VivMdvRegionReact");
        expect(config).toEqual({
            type: "SpatialDataMdvRegionReact",
            id: "c1",
            title: "r1-image",
            region: "r1",
            param: ["x", "y"],
            color_by: "cell_type",
            radius: 4,
            opacity: 0.6,
            background_filter,
            category_filters,
            contourParameter: "cell_type",
            viewState: { target: [5, 5, 0], zoom: -2 },
            defaultImageLayer: { channels: { colors: [[255, 0, 0]], contrastLimits: [[0, 100]] } },
        });
    });

    test("keeps the Viv channel settings and tone for the image layer", () => {
        const channelsStore = {
            channelsVisible: [true, false],
            colors: [
                [255, 0, 0],
                [0, 255, 0],
            ],
            contrastLimits: [
                [0, 100],
                [10, 200],
            ],
            selections: [
                { z: 0, c: 0, t: 0 },
                { z: 0, c: 2, t: 0 },
            ],
            brightness: [0.6, 0.5],
            contrast: [0.4, 0.5],
            domains: [
                [0, 255],
                [0, 255],
            ],
        };
        const config: { type: string; region: string; viv: unknown; defaultImageLayer?: unknown } = {
            type: "VivMdvRegionReact",
            region: "r1",
            viv: { channelsStore },
        };
        convertLegacyChartConfigNow(config, spatialRegion);
        expect(config.defaultImageLayer).toEqual({
            channels: {
                channelsVisible: channelsStore.channelsVisible,
                colors: channelsStore.colors,
                contrastLimits: channelsStore.contrastLimits,
                selections: channelsStore.selections,
            },
            vivLayerProps: { brightness: [0.6, 0.5], contrast: [0.4, 0.5] },
        });
    });

    test("leaves the image layer's defaults when the Viv chart saved no channels, or empty tone", () => {
        const config = { type: "VivMdvRegionReact", region: "r1", viv: { channelsStore: { brightness: [] } } };
        convertLegacyChartConfigNow(config, spatialRegion);
        expect(config).toEqual({ type: "SpatialDataMdvRegionReact", region: "r1" });
    });

    test("leaves the view to fit the data when the Viv chart saved no position", () => {
        const config = {
            type: "VivMdvRegionReact",
            region: "r1",
            viv: { channels: [] },
            viewState: { target: [1, 1] },
        };
        convertLegacyChartConfigNow(config, spatialRegion);
        expect(config).toEqual({ type: "SpatialDataMdvRegionReact", region: "r1" });
    });

    test("is converted on load only with the setting on", () => {
        const config = { type: "VivMdvRegionReact", region: "r1", viv: {} };
        expect(convertLegacyChartConfig(config, spatialRegion)).toBeUndefined();
        convertLegacyChartsSetting.set("on");
        expect(convertLegacyChartConfig(config, spatialRegion)).toBe("VivMdvRegionReact");
        expect(config).toEqual({ type: "SpatialDataMdvRegionReact", region: "r1" });
    });

    test("stays a Viv chart for a region that is not SpatialData, or that cannot be checked", () => {
        const vivOnly = withRegions({ r1: { viv_image: {} } });
        expect(
            getLegacyConversionTarget({ type: "VivMdvRegionReact", region: "r1" } as never, vivOnly),
        ).toBeUndefined();
        expect(getLegacyConversionTarget({ type: "VivMdvRegionReact", region: "r1" } as never)).toBeUndefined();
        expect(
            getLegacyConversionTarget({ type: "VivMdvRegionReact", region: "other" } as never, spatialRegion),
        ).toBeUndefined();
    });
});
