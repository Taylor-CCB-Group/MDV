import type { PickingInfo } from "@deck.gl/core";
import { ScatterplotLayer } from "@deck.gl/layers";
import { IsolatedGroupLayer } from "@spatialdata/layers";
import { describe, expect, it } from "vitest";
import { getPickingInfoWithAlternates } from "@/lib/deckPicking";
import { combineTooltipContent, getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import { deckHostLayerId } from "@/react/spatialdata/host_overlay_ids";

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
            getPointTooltip: (info) => (info ? { html: "<strong>Gene:</strong> CD4" } : null),
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

    it("finds a gate the SpatialData render stack draws, by the stack's id for it", () => {
        // the stack draws the chart's gate layer as `deck:gate_display`, not under the layer's own
        // id, and SpatialData.js's viewer draws a copy of that as `deck:gate_display-#<viewer>#`
        const gateId = deckHostLayerId("gate_display");
        const tooltipFor = (pickedLayerId: string) =>
            getCombinedScatterTooltip(pickingInfo(pickedLayerId, 0, { properties: { gateName: "g1" } }), {
                gateDisplayLayerId: gateId,
                getPointTooltip: () => null,
            });
        expect(tooltipFor(gateId)).toMatchObject({ html: expect.stringContaining("<strong>g1</strong>") });
        expect(tooltipFor(`${gateId}-#spatial-detail-x1#`)).toMatchObject({
            html: expect.stringContaining("<strong>g1</strong>"),
        });
        expect(tooltipFor(`${gateId}_other`)).toBeNull();
    });

    it("escapes text tooltip parts when combining", () => {
        expect(combineTooltipContent("A < B", { text: "C > D" })).toMatchObject({
            html: 'A &lt; B<div style="border-top: 1px solid rgba(255,255,255,0.35); margin: 6px 0;"></div>C &gt; D',
        });
    });

    it("finds the SpatialData chart's scatter inside a layer group", () => {
        const group = new IsolatedGroupLayer({
            id: "group:dev-scatter-#detail#",
            layers: [layer("deck:scatter-#detail#")],
        });
        // As the group draws it: under an id scoped to the group.
        const member = group.renderLayers().flat()[0];
        if (!(member instanceof ScatterplotLayer)) throw new Error("expected the scatter member");
        // As deck reports it once the pick has bubbled up to the group.
        const info = group.getPickingInfo({
            info: { ...pickingInfo("group:dev-scatter-#detail#", 4, 4), layer: group, sourceLayer: member },
            mode: "hover",
            sourceLayer: member,
        });

        expect(
            getCombinedScatterTooltip(info, {
                getPointTooltip: (picked) => (picked ? `row ${picked.index}` : null),
            }),
        ).toEqual({ html: "row 4" });
    });

    it("finds a gate outline drawn in a group with the points", () => {
        const group = new IsolatedGroupLayer({
            id: "shared-cell-layers",
            layers: [layer("scatter_1"), layer("gate-display")],
        });
        const gateMember = group.renderLayers().flat()[1];
        if (!(gateMember instanceof ScatterplotLayer)) throw new Error("expected the gate member");
        const info = group.getPickingInfo({
            info: {
                ...pickingInfo("shared-cell-layers", 0, { properties: { gateName: "Gate B" } }),
                layer: group,
                sourceLayer: gateMember,
            },
            mode: "hover",
            sourceLayer: gateMember,
        });

        expect(
            getCombinedScatterTooltip(info, {
                gateDisplayLayerId: "gate-display",
                getPointTooltip: () => null,
            }),
        ).toMatchObject({ html: expect.stringContaining("<strong>Gate B</strong>") });
    });
});
