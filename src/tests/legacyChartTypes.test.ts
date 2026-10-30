import BaseChart from "@/charts/BaseChart";
import "@/charts/DensityScatterPlot";
import "@/charts/WGLScatterPlot";
import "@/react/components/DeckScatterReactWrapper";
import { legacyChartsSetting } from "@/lib/browserSettings";
import { describe, expect, test } from "vitest";

const addable = (id: string) => BaseChart.types[id].allow_user_add !== false;

describe("scatter plot chart types in the Add Chart dialog", () => {
    test("legacy charts are hidden by default", () => {
        expect(legacyChartsSetting.get()).toBe("hide");
    });

    test("the deck scatter plot is the 2D scatter plot", () => {
        const type = BaseChart.types.DeckContourScatter;
        expect(type.name).toBe("2D Scatter Plot");
        expect(addable("DeckContourScatter")).toBe(true);
        expect(type.legacy).toBeFalsy();
    });

    test.each(["wgl_scatter_plot", "wgl_scatter_plot_dev", "density_scatter_plot"])("%s is legacy", (id) => {
        expect(BaseChart.types[id].legacy).toBe(true);
    });

    test("only one visible type is called 2D Scatter Plot", () => {
        const names = Object.values(BaseChart.types)
            .filter((t) => t.allow_user_add !== false && !t.legacy)
            .map((t) => t.name);
        expect(names.filter((n) => n === "2D Scatter Plot")).toHaveLength(1);
    });
});
