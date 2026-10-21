import {
    NO_COLUMN,
    adaptDotPlotConfig,
    getColorRange,
    getSplitFromDialog,
} from "@/react/components/dotplot/dotPlotConfig";
import { describe, expect, it } from "vitest";

describe("getSplitFromDialog", () => {
    const columnValues = ["ND", "R baseline", "NR baseline"];
    const none = { split_by: null, split_values: [] };

    it("keeps the chosen column and values, in the order picked", () => {
        expect(getSplitFromDialog("response_stage", ["NR baseline", "R baseline"], columnValues)).toEqual({
            split_by: "response_stage",
            split_values: ["NR baseline", "R baseline"],
        });
    });

    it("starts with no panels when no values are picked", () => {
        expect(getSplitFromDialog("response_stage", [], columnValues)).toEqual({
            split_by: "response_stage",
            split_values: [],
        });
        expect(getSplitFromDialog("response_stage", undefined, columnValues).split_values).toEqual([]);
    });

    it("drops values left over from a column chosen earlier", () => {
        expect(getSplitFromDialog("response_stage", ["Control", "R baseline"], columnValues).split_values).toEqual([
            "R baseline",
        ]);
    });

    it("means no split for None or no choice, whatever values are left over", () => {
        expect(getSplitFromDialog(NO_COLUMN, ["R baseline"], columnValues)).toEqual(none);
        expect(getSplitFromDialog("", ["R baseline"], columnValues)).toEqual(none);
        expect(getSplitFromDialog(undefined, undefined, undefined)).toEqual(none);
    });
});

const size: [number, number] = [400, 300];
const base = {
    id: "dp",
    size,
    title: "",
    legend: "",
    type: "compound_dot_plot",
    param: ["cell_type", "gene1"],
};

describe("adaptDotPlotConfig", () => {
    it("fills in defaults for a new chart", () => {
        const c = adaptDotPlotConfig(base);
        expect(c.y_axis_order).toBe("data");
        expect(c.split_by).toBeNull();
        expect(c.split_values).toEqual([]);
        expect(c.color_scale).toEqual({ palette: "reds", log: false, mode: "auto" });
        expect(c.color_legend.display).toBe(true);
        expect(c.fraction_legend.display).toBe(true);
        expect(c.axis.x).toEqual({
            size: 60,
            tickfont: 10,
            rotate_labels: true,
            fit_labels: true,
            truncate_labels: true,
            thin_labels: false,
        });
        expect(c.axis.y).toEqual({ size: 110, tickfont: 10, rotate_labels: false });
    });

    it("keeps saved settings", () => {
        const c = adaptDotPlotConfig({
            ...base,
            split_by: "condition",
            split_values: ["Disease", "Control"],
            color_scale: { palette: "blue_white_red", log: true, mode: "manual", min: 0, max: 3 },
            color_legend: { display: false },
            axis: {
                x: { size: 90, tickfont: 12, rotate_labels: false },
                y: { size: 80, tickfont: 9, rotate_labels: false },
            },
        });
        expect(c.split_by).toBe("condition");
        expect(c.split_values).toEqual(["Disease", "Control"]);
        expect(c.color_scale).toEqual({ palette: "blue_white_red", log: true, mode: "manual", min: 0, max: 3 });
        expect(c.color_legend.display).toBe(false);
        expect(c.fraction_legend.display).toBe(true);
        expect(c.axis.x.size).toBe(90);
        expect(c.axis.y.tickfont).toBe(9);
    });
});

describe("getColorRange", () => {
    const scale = { palette: "reds" as const, log: false };

    it("uses the automatic range in auto mode, even with a manual range saved", () => {
        expect(getColorRange({ ...scale, mode: "auto", min: 0, max: 5 }, [1, 2])).toEqual([1, 2]);
    });

    it("uses the manual range in manual mode", () => {
        expect(getColorRange({ ...scale, mode: "manual", min: 0, max: 5 }, [1, 2])).toEqual([0, 5]);
    });

    it("swaps a manual range entered backwards", () => {
        expect(getColorRange({ ...scale, mode: "manual", min: 5, max: 0 }, [1, 2])).toEqual([0, 5]);
    });

    it("falls back to the automatic range when the manual one is incomplete", () => {
        expect(getColorRange({ ...scale, mode: "manual", min: 0 }, [1, 2])).toEqual([1, 2]);
        expect(getColorRange({ ...scale, mode: "manual", min: 0, max: Number.NaN }, [1, 2])).toEqual([1, 2]);
    });
});
