import { describe, expect, test } from "vitest";
import {
    getDensityGridCellBounds,
    getDensityGridLayout,
    getDensityGridViewId,
    getDensityGridViewStates,
    matchesDensityGridView,
} from "./densityGridUtils";

describe("densityGridUtils", () => {
    test("lays out fields across responsive columns", () => {
        const layout = getDensityGridLayout(900, 600, 5, 260);

        expect(layout.columns).toBe(3);
        expect(layout.rows).toBe(2);
        expect(layout.cellWidth).toBe(300);
        expect(layout.cellHeight).toBe(300);
    });

    test("computes cell bounds from row and column position", () => {
        const layout = getDensityGridLayout(900, 600, 5, 260);

        expect(getDensityGridCellBounds(layout, 4)).toEqual({
            x: 300,
            y: 300,
            width: 300,
            height: 300,
        });
    });

    test("creates stable view ids and scopes sublayers", () => {
        const viewId = getDensityGridViewId("chart:#1", "marker/a", 2);

        expect(viewId).toBe("density-grid-chart__1-2-marker_a");
        expect(matchesDensityGridView(`${viewId}-weights`, viewId)).toBe(true);
        expect(matchesDensityGridView(`${viewId}-weights`, "other-view")).toBe(false);
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
