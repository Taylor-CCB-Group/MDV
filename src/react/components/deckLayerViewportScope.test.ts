import { describe, expect, test } from "vitest";
import { getDensityGridViewId } from "./densityGridUtils";
import {
    MDV_DECK_LAYER_VIEWPORT_SCOPE,
    shouldDrawDeckLayerInViewport,
} from "./deckLayerViewportScope";

describe("deckLayerViewportScope", () => {
    const detailId = "chart-1detail-react";
    const vivSuffix = `-#${detailId}#`;
    const gridView = getDensityGridViewId("chart-1", "field_a", 0);

    test("chart-shared layers draw on overlay detail and every grid cell", () => {
        const gate = {
            id: `gate_${vivSuffix}`,
            props: { [MDV_DECK_LAYER_VIEWPORT_SCOPE]: "chart-shared" },
        };
        const json = {
            id: `json_${vivSuffix}`,
            props: { [MDV_DECK_LAYER_VIEWPORT_SCOPE]: "chart-shared" },
        };

        expect(shouldDrawDeckLayerInViewport(gate, detailId, { overlayDetailVivId: vivSuffix })).toBe(true);
        expect(shouldDrawDeckLayerInViewport(json, detailId, { overlayDetailVivId: vivSuffix })).toBe(true);
        expect(shouldDrawDeckLayerInViewport(gate, gridView, { overlayDetailVivId: vivSuffix })).toBe(true);
        expect(shouldDrawDeckLayerInViewport(json, gridView, { overlayDetailVivId: vivSuffix })).toBe(true);
    });

    test("per-viewport layers only draw in their grid cell", () => {
        const density = {
            id: `${gridView}-density`,
            props: {
                [MDV_DECK_LAYER_VIEWPORT_SCOPE]: "per-viewport",
                viewId: gridView,
            },
        };
        const otherCell = getDensityGridViewId("chart-1", "field_b", 1);

        expect(shouldDrawDeckLayerInViewport(density, gridView, { overlayDetailVivId: vivSuffix })).toBe(true);
        expect(shouldDrawDeckLayerInViewport(density, otherCell, { overlayDetailVivId: vivSuffix })).toBe(false);
        expect(shouldDrawDeckLayerInViewport(density, detailId, { overlayDetailVivId: vivSuffix })).toBe(false);
    });

    test("infers chart-shared for selection layers without explicit scope", () => {
        const selection = { id: "selection_-#chartdetail-react#", props: {} };
        expect(shouldDrawDeckLayerInViewport(selection, gridView)).toBe(true);
    });
});
