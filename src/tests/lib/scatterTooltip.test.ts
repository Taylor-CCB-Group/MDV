import { ScatterplotLayer } from "@deck.gl/layers";
import type { PickingInfo } from "@deck.gl/core";
import { describe, expect, it } from "vitest";
import { getPickingInfoWithAlternates } from "@/lib/deckPicking";
import { combineTooltipContent, getCombinedScatterTooltip } from "@/lib/scatterTooltip";

function layer(id: string) {
    return new ScatterplotLayer({ id, data: [] });
}

function pickingInfo(id: string, index: number, object?: unknown): PickingInfo {
    return {
        color: null,
        index,
        layer: layer(id),
        object,
        picked: true,
        pixelRatio: 1,
        x: 12,
        y: 34,
    };
}

describe("scatterTooltip", () => {
    it("combines gate context with point details", () => {
        const pointInfo = pickingInfo("scatter_1", 3, 21);
        const gateInfo = pickingInfo("gate-display", 1, {
            properties: { gateName: "Gate <A>" },
        });
        const richInfo = getPickingInfoWithAlternates(pointInfo, {
            pickMultipleObjects: () => [gateInfo],
        });

        const tooltip = getCombinedScatterTooltip(richInfo, {
            gateDisplayLayerId: "gate-display",
            getPointTooltip: (info) => info ? { html: "<strong>Gene:</strong> CD4" } : null,
        });

        expect(tooltip).toMatchObject({
            html: expect.stringContaining("<strong>Gate &lt;A&gt;</strong>"),
        });
        expect(tooltip).toMatchObject({
            html: expect.stringContaining("Click on the label to edit"),
        });
        expect(tooltip).toMatchObject({
            html: expect.stringContaining("border-top"),
        });
        expect(tooltip).toMatchObject({
            html: expect.stringContaining("<strong>Gene:</strong> CD4"),
        });
    });

    it("escapes text tooltip parts when combining", () => {
        expect(combineTooltipContent("A < B", { text: "C > D" })).toMatchObject({
            html: "A &lt; B<div style=\"border-top: 1px solid rgba(255,255,255,0.35); margin: 6px 0;\"></div>C &gt; D",
        });
    });
});
