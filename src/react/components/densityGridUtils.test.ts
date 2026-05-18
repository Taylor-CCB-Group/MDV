import { describe, expect, test } from "vitest";
import { getVivGridDetailViewId } from "./chartArrayGridUtils";
import {
    getDensityGridViewId,
    getDensityGridViewStates,
    getViewportAtCanvasPoint,
    hasUsableOrthographicViewState,
    isDensityGridViewport,
    isEditableSelectionLayerId,
    matchesDensityGridView,
    shouldDrawLayerInDeckDensityGrid,
    shouldDrawLayerInViewport,
    supportsDensityGridMode,
    unprojectCanvasPoint,
} from "./densityGridUtils";

describe("densityGridUtils", () => {
    test("creates stable view ids and scopes sublayers", () => {
        const viewId = getDensityGridViewId("chart:#1", "marker/a", 2);

        expect(viewId).toBe("density-grid-chart__1-2-marker_a");
        expect(matchesDensityGridView(`${viewId}-weights`, viewId)).toBe(true);
        expect(matchesDensityGridView(`${viewId}-selection`, viewId)).toBe(true);
        expect(matchesDensityGridView(`${viewId}-weights`, "other-view")).toBe(false);
    });

    test("supports density grid for deck and viv spatial chart types", () => {
        expect(supportsDensityGridMode("DeckContourScatter")).toBe(true);
        expect(supportsDensityGridMode("DeckDensity")).toBe(true);
        expect(supportsDensityGridMode("VivMdvRegionReact")).toBe(true);
        expect(supportsDensityGridMode("viv_scatter_plot")).toBe(false);
        expect(supportsDensityGridMode("wgl_scatter_plot")).toBe(false);
        expect(supportsDensityGridMode(undefined)).toBe(false);
    });

    test("creates viv detail view ids from grid view ids", () => {
        const detailId = getVivGridDetailViewId("chart:#1", "marker/a", 2);
        expect(detailId).toBe("density-grid-chart__1-2-marker_adetail-react");
        expect(matchesDensityGridView(`${detailId}-gate`, detailId)).toBe(true);
    });

    test("detects usable orthographic view state", () => {
        expect(hasUsableOrthographicViewState(undefined)).toBe(false);
        expect(hasUsableOrthographicViewState({ target: [1, 2, 0], zoom: 3 })).toBe(true);
        expect(hasUsableOrthographicViewState({ target: [Number.NaN, 2, 0], zoom: 3 })).toBe(false);
    });

    test("identifies density grid viewports", () => {
        expect(isDensityGridViewport("density-grid-chart__1-0-field")).toBe(true);
        expect(isDensityGridViewport("my-chartdetail-react")).toBe(false);
    });

    test("draws overlay gates on detail viewports after density grid mode", () => {
        const detailId = "chart-1detail-react";
        const vivSuffix = `-#${detailId}#`;
        const gateLayer = { id: `gate_${vivSuffix}`, props: {} };
        const staleGridGate = {
            id: `${getDensityGridViewId("chart-1", "field_a", 0)}-gate`,
            props: { viewId: getDensityGridViewId("chart-1", "field_a", 0) },
        };

        expect(shouldDrawLayerInViewport(gateLayer, detailId, vivSuffix)).toBe(true);
        expect(shouldDrawLayerInViewport(staleGridGate, detailId, vivSuffix)).toBe(false);
        expect(
            shouldDrawLayerInViewport(
                staleGridGate,
                getDensityGridViewId("chart-1", "field_a", 0),
                `-#${getDensityGridViewId("chart-1", "field_a", 0)}#`,
            ),
        ).toBe(true);
    });

    test("identifies editable selection layers", () => {
        expect(isEditableSelectionLayerId("selection_-#chartdetail-react#")).toBe(true);
        expect(isEditableSelectionLayerId("gate_-#chartdetail-react#")).toBe(false);
    });

    test("draws selection in every density grid viewport", () => {
        const selection = { id: "selection_-#chartdetail-react#", props: {} };
        const gridView = getDensityGridViewId("chart-1", "field_a", 0);

        expect(shouldDrawLayerInDeckDensityGrid(selection, gridView)).toBe(true);
        expect(shouldDrawLayerInViewport(selection, gridView, "-#ignored#")).toBe(true);
        expect(shouldDrawLayerInDeckDensityGrid(selection, "chart-1detail-react")).toBe(false);
    });

    test("unprojects through the sub-viewport under the pointer", () => {
        const viewports = [
            {
                id: "density-grid-a-0-x",
                x: 0,
                y: 0,
                width: 100,
                height: 100,
                unproject: ([px, py]: number[]) => [px, py, 0],
            },
            {
                id: "density-grid-a-1-y",
                x: 100,
                y: 0,
                width: 100,
                height: 100,
                unproject: ([px, py]: number[]) => [px + 1000, py, 0],
            },
        ];

        expect(getViewportAtCanvasPoint(viewports, 10, 10)?.id).toBe("density-grid-a-0-x");
        expect(getViewportAtCanvasPoint(viewports, 150, 10)?.id).toBe("density-grid-a-1-y");
        expect(unprojectCanvasPoint(viewports, 150, 20)).toEqual([1050, 20, 0]);
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
