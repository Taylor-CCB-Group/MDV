import { describe, expect, test } from "vitest";
import {
    adaptSplatterConfig,
    getSharedSplatterViewState,
    getSplatterLayout,
    getVisibleSplatterCategories,
    type SplatterPlotConfig,
} from "./splatterPlotUtils";

describe("SplatterPlotReactWrapper helpers", () => {
    test("adaptSplatterConfig separates density fields from coordinate params", () => {
        const config = adaptSplatterConfig({
            id: "splatter-1",
            size: [300, 200],
            title: "splatter",
            legend: "",
            type: "DeckSplatter",
            param: ["x", "y", "field_a", "field_b"],
            category: "cell_type",
        } as SplatterPlotConfig);

        expect(config.param).toEqual(["x", "y"]);
        expect(config.densityFields).toEqual(["field_a", "field_b"]);
        expect(config.contour_fill).toBe(true);
    });

    test("getVisibleSplatterCategories keeps category order and drops empty rows", () => {
        const categories = getVisibleSplatterCategories(
            {
                field: "cell_type",
                name: "Cell Type",
                datatype: "text",
                values: ["T cell", "B cell", "NK"],
                data: new Uint8Array([0, 1, 0, 2, 2]),
            } as any,
            new Uint32Array([0, 2, 4]),
        );

        expect(categories.map((category) => category.label)).toEqual(["T cell", "NK"]);
        expect(Array.from(categories[0].rows)).toEqual([0, 2]);
        expect(Array.from(categories[1].rows)).toEqual([4]);
    });

    test("getSplatterLayout reserves label/header space and splits remaining plot evenly", () => {
        const layout = getSplatterLayout(640, 480, ["A", "much longer category"], ["Field A", "Field B"]);

        expect(layout.labelWidth).toBeGreaterThan(88);
        expect(layout.headerHeight).toBeGreaterThanOrEqual(58);
        expect(layout.cellWidth).toBeCloseTo(layout.plotWidth / 2);
        expect(layout.cellHeight).toBeCloseTo(layout.plotHeight / 2);
    });

    test("getSharedSplatterViewState fits the filtered rows inside each cell", () => {
        const viewState = getSharedSplatterViewState(
            {
                data: new Float32Array([0, 10, 20]),
            } as any,
            {
                data: new Float32Array([0, 5, 10]),
            } as any,
            new Uint32Array([0, 2]),
            120,
            80,
        );

        expect(viewState.target).toEqual([10, 5, 0]);
        expect(viewState.zoom).toBeGreaterThan(1);
        expect(viewState.minZoom).toBe(-50);
    });
});
