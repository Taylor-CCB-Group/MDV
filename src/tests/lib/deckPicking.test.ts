import { ScatterplotLayer } from "@deck.gl/layers";
import type { PickingInfo } from "@deck.gl/core";
import { describe, expect, it, vi } from "vitest";
import {
    findPickedInfo,
    getPickedInfos,
    getPickingInfoWithAlternates,
    pickingInfoMatchesLayer,
} from "@/lib/deckPicking";

function layer(id: string) {
    return new ScatterplotLayer({ id, data: [] });
}

function pickingInfo(id: string, index: number, object?: unknown, sourceId?: string): PickingInfo {
    return {
        color: null,
        index,
        layer: layer(id),
        object,
        picked: true,
        pixelRatio: 1,
        sourceLayer: sourceId ? layer(sourceId) : null,
        x: 12,
        y: 34,
    };
}

describe("deckPicking", () => {
    it("wraps the primary picking info when no deck picker is available", () => {
        const primaryInfo = pickingInfo("scatter_1", 2, 42);

        const richInfo = getPickingInfoWithAlternates(primaryInfo, null);

        expect(getPickedInfos(richInfo)).toEqual([primaryInfo]);
    });

    it("accumulates alternate picks and deduplicates the primary hit", () => {
        const primaryInfo = pickingInfo("scatter_1", 2, 42);
        const duplicatePrimaryInfo = pickingInfo("scatter_1", 2, 42);
        const gateInfo = pickingInfo("gate-display", 1, { properties: { gateName: "A" } });
        const pickMultipleObjects = vi.fn(() => [duplicatePrimaryInfo, gateInfo]);

        const richInfo = getPickingInfoWithAlternates(
            primaryInfo,
            { pickMultipleObjects },
            { depth: 4, radius: 3 },
        );

        expect(pickMultipleObjects).toHaveBeenCalledWith({ x: 12, y: 34, radius: 3, depth: 4 });
        expect(getPickedInfos(richInfo)).toEqual([primaryInfo, gateInfo]);
    });

    it("finds alternate picking info by layer predicate", () => {
        const primaryInfo = pickingInfo("scatter_1", 2, 42);
        const gateInfo = pickingInfo("gate-display", 1, { properties: { gateName: "A" } });
        const richInfo = getPickingInfoWithAlternates(primaryInfo, {
            pickMultipleObjects: () => [gateInfo],
        });

        expect(findPickedInfo(richInfo, (info) => info.layer?.id === "gate-display")).toBe(gateInfo);
    });

    it("matches source layer ids for composite layer picks", () => {
        const info = pickingInfo("spatial", 3, 84, "spatial.scatterplot-points");

        expect(
            pickingInfoMatchesLayer(info, (layerId) => layerId.includes("spatial.scatterplot")),
        ).toBe(true);
    });
});
