import type { PickingInfo } from "@deck.gl/core";
import { ScatterplotLayer } from "@deck.gl/layers";
import type { MjolnirKeyEvent } from "mjolnir.js";
import { describe, expect, it, vi } from "vitest";
import { SHARED_CELL_LAYERS_ID, sharedCellLayersGroup } from "./densityGridSharedLayers";
import { shouldDrawInDensityGridCell } from "./densityGridUtils";

/** A pick of `member` as deck hands it to the group: bubbled up, naming the group. */
function pickThrough(group: ReturnType<typeof sharedCellLayersGroup>, memberIndex: number): PickingInfo {
    const member = group.renderLayers().flat()[memberIndex];
    if (!(member instanceof ScatterplotLayer)) throw new Error("expected a scatter member");
    return group.getPickingInfo({
        info: { color: null, index: 2, layer: group, sourceLayer: member, picked: true, pixelRatio: 1, x: 1, y: 1 },
        mode: "hover",
        sourceLayer: member,
    });
}

/** Any event will do: the handlers only pass it on. */
function anEvent(): MjolnirKeyEvent {
    const element = document.createElement("div");
    return {
        type: "keyup",
        key: "x",
        srcEvent: new KeyboardEvent("keyup"),
        target: element,
        rootElement: element,
        handled: false,
        stopPropagation: () => {},
        stopImmediatePropagation: () => {},
    };
}

describe("sharedCellLayersGroup", () => {
    it("draws in every cell", () => {
        expect(shouldDrawInDensityGridCell(SHARED_CELL_LAYERS_ID, "density-grid-chart-field-0")).toBe(true);
    });

    it("passes a pick's events on to the member it came through", () => {
        const onClick = vi.fn(() => true);
        const onHover = vi.fn();
        const otherClick = vi.fn();
        const group = sharedCellLayersGroup([
            new ScatterplotLayer({ id: "other", data: [], onClick: otherClick }),
            new ScatterplotLayer({ id: "points", data: [], onClick, onHover }),
        ]);
        const info = pickThrough(group, 1);
        const event = anEvent();

        expect(group.props.onClick?.(info, event)).toBe(true);
        expect(onClick).toHaveBeenCalledWith(info, event);
        expect(otherClick).not.toHaveBeenCalled();
        group.props.onHover?.(info, event);
        expect(onHover).toHaveBeenCalledWith(info, event);
    });

    it("leaves a pick unhandled when its member has no handler", () => {
        const group = sharedCellLayersGroup([new ScatterplotLayer({ id: "points", data: [] })]);
        expect(group.props.onDragStart?.(pickThrough(group, 0), anEvent())).toBe(false);
    });
});
