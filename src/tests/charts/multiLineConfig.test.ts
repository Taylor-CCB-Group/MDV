import { describe, expect, it } from "vitest";
import { MultiLineChartConfigSchema, MultiLineReactConfigSchema } from "@/charts/schemas/ChartConfigSchema";
import {
    adaptMultiLineConfig,
    bandWidthFromFraction,
    getBandWidthFraction,
    getBandwidth,
    getBrushFilterRange,
    getCellXDomains,
    getXDomain,
    keepBrushes,
    type SavedMultiLineConfig,
    withBrush,
} from "@/react/components/multiline/multiLineConfig";

/** One of a marker QC view's multi-line charts, as the legacy chart saved it. */
const legacy = {
    id: "a",
    size: [300, 300] as [number, number],
    title: "DAPI",
    legend: "",
    type: "multi_line_chart",
    param: ["zscore_clipped_DAPI", "sample_id"],
    stacked: false,
    fill: false,
    scaletrim: "0.05",
    band_width: 0.1,
    intervals: 100,
    axis: { x: { size: 30, label: "DAPI", textsize: 13 }, y: { size: 45, label: "density", textsize: 13 } },
    color_legend: { display: true, pos: [10, 20] as [number, number] },
} satisfies SavedMultiLineConfig;

describe("adaptMultiLineConfig", () => {
    it("moves a legacy param to categories and values", () => {
        const config = adaptMultiLineConfig(legacy);
        expect(config.param).toEqual([]);
        expect(config.values).toEqual(["zscore_clipped_DAPI"]);
        expect(config.categories).toBe("sample_id");
        expect(config.scaletrim).toBe("0.05");
        expect(config.band_width).toBe(0.1);
        expect(config.intervals).toBe(100);
        expect(config.color_legend).toEqual({ display: true, pos: [10, 20] });
        expect(config.axis.x).toEqual({ size: 30, tickfont: 13, rotate_labels: false, label: "DAPI" });
    });

    it("keeps named props over param", () => {
        const config = adaptMultiLineConfig({ ...legacy, param: [], categories: "batch", values: ["CD3"] });
        expect(config.categories).toBe("batch");
        expect(config.values).toEqual(["CD3"]);
    });

    it("reads 'none' or no trim as untrimmed, and hides the legend when stacked by default", () => {
        const { color_legend: _, ...rest } = legacy;
        const config = adaptMultiLineConfig({ ...rest, scaletrim: "none", stacked: true });
        expect(config.scaletrim).toBeNull();
        expect(config.color_legend.display).toBe(false);
    });

    it("gives a config both schemas accept", () => {
        expect(MultiLineChartConfigSchema.safeParse(legacy).success).toBe(true);
        const config = adaptMultiLineConfig(legacy);
        expect(MultiLineReactConfigSchema.safeParse({ ...config, type: "multi_line_chart_react" }).success).toBe(true);
    });

    it("accepts a legacy config saved with a null trim, and reads it as untrimmed", () => {
        const untrimmed = { ...legacy, scaletrim: null };
        expect(MultiLineChartConfigSchema.safeParse(untrimmed).success).toBe(true);
        expect(adaptMultiLineConfig(untrimmed).scaletrim).toBeNull();
    });
});

describe("several values", () => {
    it("keeps several values, and reads only a known grid", () => {
        const config = adaptMultiLineConfig({ ...legacy, param: [], categories: "sample_id", values: ["CD3", "CD4"] });
        expect(config.values).toEqual(["CD3", "CD4"]);
        expect(config.grid).toBeNull();
        expect(config.share_x).toBe(false);
        expect(adaptMultiLineConfig({ ...legacy, grid: { by: "values" } }).grid).toEqual({ by: "values" });
        expect(adaptMultiLineConfig({ ...legacy, grid: { by: "densityFields" } }).grid).toBeNull();
    });

    it("gives a config with brushes and a grid the schema accepts", () => {
        const config = adaptMultiLineConfig({
            ...legacy,
            values: ["CD3", "CD4"],
            grid: { by: "values" },
            share_x: true,
            brush: { CD3: [0, 1] },
        });
        expect(MultiLineReactConfigSchema.safeParse({ ...config, type: "multi_line_chart_react" }).success).toBe(true);
    });
});

describe("cell x domains", () => {
    const a = { minMax: [-2, 8] as [number, number], quantiles: { "0.05": [0, 4] as [number, number] } };
    const b = { minMax: [-5, 3] as [number, number], quantiles: { "0.05": [-1, 2] as [number, number] } };

    it("gives each cell its own domain by default", () => {
        expect(getCellXDomains([a, b], "0.05", false)).toEqual([
            { domain: [0, 4], clamp: true },
            { domain: [-1, 2], clamp: true },
        ]);
    });

    it("shares the union of the domains when asked", () => {
        const [first, second] = getCellXDomains([a, b], "0.05", true);
        expect(first).toEqual({ domain: [-1, 4], clamp: true });
        // the same object, so cells key their densities by it alike
        expect(second).toBe(first);
        expect(getCellXDomains([a, b], null, true)[0]).toEqual({ domain: [-5, 8], clamp: false });
    });

    it("has nothing to share with one value", () => {
        expect(getCellXDomains([a], "0.05", true)).toEqual([{ domain: [0, 4], clamp: true }]);
    });
});

describe("x domain and bandwidth", () => {
    const range = { minMax: [-2, 8] as [number, number], quantiles: { "0.05": [0, 4] as [number, number] } };

    it("trims to the quantile range and clamps", () => {
        expect(getXDomain(range, "0.05")).toEqual({ domain: [0, 4], clamp: true });
        expect(getXDomain(range, null)).toEqual({ domain: [-2, 8], clamp: false });
        // no quantiles for that trim: the full range
        expect(getXDomain(range, "0.01")).toEqual({ domain: [-2, 8], clamp: false });
    });

    it("reads band_width as a fraction of the trimmed range, and in data units otherwise", () => {
        expect(getBandwidth({ band_width: 0.1 }, [0, 4], true)).toBeCloseTo(0.4);
        expect(getBandwidth({ band_width: null }, [0, 4], true)).toBeCloseTo(0.4);
        expect(getBandwidth({ band_width: 0.5 }, [-2, 8], false)).toBe(0.5);
        expect(getBandwidth({ band_width: null }, [-2, 8], false)).toBeCloseTo(0.1);
    });

    it("shows band_width as a fraction of the range on the slider, and converts it back", () => {
        expect(getBandWidthFraction({ band_width: 0.2 }, [0, 4], true)).toBe(0.2);
        expect(getBandWidthFraction({ band_width: 0.5 }, [-2, 8], false)).toBeCloseTo(0.05);
        expect(getBandWidthFraction({ band_width: null }, [-2, 8], false)).toBeCloseTo(0.01);
        expect(bandWidthFromFraction(0.2, [0, 4], true)).toBe(0.2);
        expect(bandWidthFromFraction(0.05, [-2, 8], false)).toBeCloseTo(0.5);
        // a fraction gives the same kernel trimmed or not
        expect(getBandwidth({ band_width: bandWidthFromFraction(0.1, [0, 4], true) }, [0, 4], true)).toBeCloseTo(
            getBandwidth({ band_width: bandWidthFromFraction(0.1, [0, 4], false) }, [0, 4], false),
        );
    });
});

describe("brush", () => {
    it("is read by field from a saved config, and left empty otherwise", () => {
        expect(adaptMultiLineConfig(legacy).brush).toEqual({});
        expect(adaptMultiLineConfig({ ...legacy, brush: { CD3: [0.2, 1.5], CD4: [1] } }).brush).toEqual({
            CD3: [0.2, 1.5],
        });
    });

    it("reads a single saved range as the first value's", () => {
        expect(adaptMultiLineConfig({ ...legacy, brush: [0.2, 1.5] }).brush).toEqual({
            zscore_clipped_DAPI: [0.2, 1.5],
        });
        expect(adaptMultiLineConfig({ ...legacy, brush: [0.2] as never }).brush).toEqual({});
    });

    it("is set and cleared per field", () => {
        const one = withBrush({}, "CD3", [0, 1]);
        const two = withBrush(one, "CD4", [2, 3]);
        expect(two).toEqual({ CD3: [0, 1], CD4: [2, 3] });
        expect(withBrush(two, "CD3", null)).toEqual({ CD4: [2, 3] });
        expect(keepBrushes(two, ["CD4", "CD8"])).toEqual({ CD4: [2, 3] });
    });

    it("filters the brushed range, in either direction", () => {
        expect(getBrushFilterRange([0.5, 1.5], [0, 4])).toEqual({ min: 0.5, max: 1.5 });
        expect(getBrushFilterRange([1.5, 0.5], [0, 4])).toEqual({ min: 0.5, max: 1.5 });
    });

    it("takes in everything beyond an end of the axis that the brush reaches", () => {
        expect(getBrushFilterRange([0, 1.5], [0, 4])).toEqual({ min: Number.NEGATIVE_INFINITY, max: 1.5 });
        expect(getBrushFilterRange([2, 4], [0, 4])).toEqual({ min: 2, max: Number.POSITIVE_INFINITY });
    });
});
