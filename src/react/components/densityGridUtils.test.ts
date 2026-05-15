import { describe, expect, test } from "vitest";
import {
    getDensityGridCellBounds,
    getDensityGridLayout,
    getDensityGridViewId,
    getDensityGridViewStates,
    getDensityGridVisibleCellIndices,
    matchesDensityGridView,
    supportsDensityGridMode,
} from "./densityGridUtils";

describe("densityGridUtils", () => {
    test("lays out fields across responsive columns with square cells", () => {
        const layout = getDensityGridLayout(900, 600, 5, 260);

        expect(layout.columns).toBe(3);
        expect(layout.rows).toBe(2);
        expect(layout.cellSize).toBe(260);
        expect(layout.cellWidth).toBe(260);
        expect(layout.cellHeight).toBe(260);
        expect(layout.contentWidth).toBe(780);
        expect(layout.contentHeight).toBe(520);
    });

    test("content can exceed viewport height for scrolling", () => {
        const layout = getDensityGridLayout(900, 400, 5, 260);

        expect(layout.contentHeight).toBe(520);
        expect(layout.contentHeight).toBeGreaterThan(400);
    });

    test("collects visible cell indices from virtualized row and column ranges", () => {
        const layout = getDensityGridLayout(900, 600, 5, 260);

        expect(
            getDensityGridVisibleCellIndices(
                layout,
                5,
                [{ index: 0 }],
                [{ index: 0 }, { index: 1 }],
            ),
        ).toEqual([0, 1]);
        expect(
            getDensityGridVisibleCellIndices(
                layout,
                5,
                [{ index: 1 }],
                [{ index: 0 }, { index: 1 }],
            ),
        ).toEqual([3, 4]);
        expect(
            getDensityGridVisibleCellIndices(
                layout,
                4,
                [{ index: 1 }],
                [{ index: 0 }, { index: 1 }, { index: 2 }],
            ),
        ).toEqual([3]);
    });

    test("computes cell bounds from row and column position", () => {
        const layout = getDensityGridLayout(900, 600, 5, 260);

        expect(getDensityGridCellBounds(layout, 4)).toEqual({
            x: 260,
            y: 260,
            width: 260,
            height: 260,
        });
    });

    test("creates stable view ids and scopes sublayers", () => {
        const viewId = getDensityGridViewId("chart:#1", "marker/a", 2);

        expect(viewId).toBe("density-grid-chart__1-2-marker_a");
        expect(matchesDensityGridView(`${viewId}-weights`, viewId)).toBe(true);
        expect(matchesDensityGridView(`${viewId}-selection`, viewId)).toBe(true);
        expect(matchesDensityGridView(`${viewId}-weights`, "other-view")).toBe(false);
    });

    test("supports density grid only for contour deck chart types", () => {
        expect(supportsDensityGridMode("DeckContourScatter")).toBe(true);
        expect(supportsDensityGridMode("DeckDensity")).toBe(true);
        expect(supportsDensityGridMode("wgl_scatter_plot")).toBe(false);
        expect(supportsDensityGridMode(undefined)).toBe(false);
    });

    test("maps shared view state onto each view id", () => {
        const viewStates = getDensityGridViewStates(
            ["view-a", "view-b"],
            { target: [1, 2, 0], zoom: 3, minZoom: -50 },
        );

        expect(viewStates["view-a"]).toMatchObject({ target: [1, 2, 0], zoom: 3 });
        expect(viewStates["view-b"]).toMatchObject({ target: [1, 2, 0], zoom: 3 });
    });
});
