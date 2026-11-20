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
import { DEFAULT_ROTATION_ORBIT, DEFAULT_ROTATION_X } from "@/react/orbitView";

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
            for (const type of ["VivMdvRegionReact", "DeckContourScatter", "DeckScatter3D"]) {
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

describe("the classic 3D scatter plot", () => {
    // saved by the classic chart (from real views)
    const classic = () => ({
        type: "wgl_3d_scatter_plot",
        id: "c3",
        title: "PCA",
        param: ["pc1", "pc2", "pc3"],
        camera: { distance: 63.96, theta: -1.038, phi: 0.261 },
        center: [0.13, 0.19, 0.92],
        axis_scales: [1, 1, 1],
        radius: 5,
        opacity: 0.8,
        size: [300, 300],
        brush: "default",
        on_filter: "hide",
        default_color: "#377eb8",
        color_by: "leiden",
        tooltip: { show: false },
    });

    /** Column ranges, as `getMinMaxForColumn` gives them. */
    const withRanges = (ranges: Record<string, [number, number]>) =>
        ({ getMinMaxForColumn: (column: string) => ranges[column] }) as unknown as DataStore;

    test("loads as the deck 3D scatter plot whatever the setting", () => {
        convertLegacyChartsSetting.set("off");
        const config: Record<string, unknown> & { type: string } = classic();
        expect(convertLegacyChartConfig(config)).toBe("wgl_3d_scatter_plot");
        expect(config).toMatchObject({
            type: "DeckScatter3D",
            dimension: "3d",
            orbit_axis: "Y",
            param: ["pc1", "pc2", "pc3"],
            title: "PCA",
            opacity: 0.8,
            on_filter: "hide",
            color_by: "leiden",
            tooltip: { show: false },
            axis_scales: [1, 1, 1],
        });
        for (const key of ["camera", "center", "brush"]) expect(config).not.toHaveProperty(key);
    });

    test("keeps its camera, as an orbit view with y up", () => {
        const config: { type: string; viewState?: any } = classic();
        convertLegacyChartConfigNow(config);
        const { target, zoom, rotationOrbit, rotationX } = config.viewState;
        expect(target).toEqual([0.13, 0.19, 0.92]);
        // a quarter turn on from theta, and phi up
        expect(rotationOrbit).toBeCloseTo((-1.038 * 180) / Math.PI + 90);
        expect(rotationX).toBeCloseTo((0.261 * 180) / Math.PI);
        // a data unit is (plot height) / (distance * tan(22.5°)) pixels, as the classic chart drew it
        expect(2 ** zoom).toBeCloseTo(250 / (63.96 * Math.tan(Math.PI / 8)), 1);
    });

    test("the classic default view is the one a new chart starts from", () => {
        const config: { type: string; viewState?: any } = classic();
        convertLegacyChartConfigNow(config);
        expect(config.viewState.rotationOrbit).toBeCloseTo(DEFAULT_ROTATION_ORBIT, 0);
        expect(config.viewState.rotationX).toBeCloseTo(DEFAULT_ROTATION_X, 0);
    });

    test("replaces a view saved when this type was briefly the deck chart", () => {
        const config: ReturnType<typeof classic> & { viewState?: any } = {
            ...classic(),
            camera: { distance: 49.6, phi: 0, theta: 0 },
            center: [0, 0, 0],
            viewState: { target: [8.3, 2.6, 5.7], zoom: 3.8, rotationOrbit: 0, rotationX: 0 },
        };
        convertLegacyChartConfigNow(config);
        expect(config.viewState.target).toEqual([0, 0, 0]);
        expect(config.viewState.rotationOrbit).toBeCloseTo(90);
    });

    test("is fitted to its data when it saved no camera", () => {
        const { camera, center, ...config }: Record<string, unknown> & { type: string } = {
            ...classic(),
            viewState: { target: [1, 2, 3], zoom: 1 },
        };
        convertLegacyChartConfigNow(config);
        expect(config).not.toHaveProperty("viewState");
    });

    test("keeps its axis scales", () => {
        const config = { ...classic(), axis_scales: [0.5, 0.5, 3] };
        convertLegacyChartConfigNow(config);
        expect(config.axis_scales).toEqual([0.5, 0.5, 3]);
    });

    test("draws its points the size they were at the saved view", () => {
        const config: { type: string; radius: number; course_radius?: number; viewState?: any } = classic();
        const ranges = { pc1: [-10, 30], pc2: [-5, 15], pc3: [-8, 8] } satisfies Record<string, [number, number]>;
        convertLegacyChartConfigNow(config, withRanges(ranges));
        // the deck chart's radius is a fraction of the larger x / y range, in data units
        const deckPixelRadius =
            ((config.radius * (config.course_radius ?? 1)) / 10000) * 40 * 2 ** config.viewState.zoom;
        // the classic radius is a diameter
        expect(deckPixelRadius).toBeCloseTo(2.5, 0);
    });

    test("does not start hiding rows with filters the classic chart never applied", () => {
        const config: Record<string, unknown> & { type: string } = {
            ...classic(),
            background_filter: { column: "cell_type", categories: ["T", "B"], category: "T" },
            category_filters: [{ column: "sample", category: ["a"] }],
        };
        convertLegacyChartConfigNow(config);
        expect(config).not.toHaveProperty("background_filter");
        expect(config).not.toHaveProperty("category_filters");
    });

    test("sizes points for data that only varies along z", () => {
        const config: { type: string; radius: number; course_radius?: number; viewState?: any } = classic();
        convertLegacyChartConfigNow(config, withRanges({ pc1: [2, 2], pc2: [3, 3], pc3: [-8, 8] }));
        const deckPixelRadius =
            ((config.radius * (config.course_radius ?? 1)) / 10000) * 16 * 2 ** config.viewState.zoom;
        expect(deckPixelRadius).toBeCloseTo(2.5, 0);
    });

    test("keeps its point size when the column ranges are not known", () => {
        const config = classic();
        convertLegacyChartConfigNow(config);
        expect(config.radius).toBe(5);
    });

    test("the deck chart under its older name converts too", () => {
        const config = { type: "wgl_3d_scatter_plot_dev", viewState: { target: [0, 0, 0], zoom: 1 } };
        expect(convertLegacyChartConfig(config)).toBe("wgl_3d_scatter_plot_dev");
        expect(config).toEqual({ type: "DeckScatter3D", dimension: "3d", viewState: { target: [0, 0, 0], zoom: 1 } });
    });
});
