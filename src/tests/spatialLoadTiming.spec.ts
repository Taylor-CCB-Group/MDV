import { renderHook } from "@testing-library/react";
import type { DeckGLRef } from "deck.gl";
import { beforeEach, describe, expect, it, vi } from "vitest";

const recordSpatialLoaded = vi.fn();
vi.mock("@/react/spatialdata/perf", () => ({ recordSpatialLoaded }));

const { useMarkSpatialChartOpened, useSpatialLoadTiming } = await import("@/react/spatialdata/load_timing");

/** A deck at `zoom` whose layers report loaded or not, as `layer.isLoaded` does. */
function deckAt(zoom: number, ...loaded: boolean[]) {
    const viewport = { width: 800, height: 600, zoom, center: [0, 0, 0] };
    const deck = { props: { layers: [loaded.map((isLoaded) => ({ isLoaded }))] }, getViewports: () => [viewport] };
    return { deck } as unknown as DeckGLRef;
}

const describeLoad = () => ({ layers: "1", image_layers: "1" });

describe("useSpatialLoadTiming", () => {
    let now: number;
    let deckRef: { current: DeckGLRef | null };
    let chart: object;
    let afterRender: () => void;

    beforeEach(() => {
        recordSpatialLoaded.mockClear();
        now = 1000;
        vi.spyOn(performance, "now").mockImplementation(() => now);
        deckRef = { current: null };
        chart = {};
        renderHook(() => useMarkSpatialChartOpened(chart));
        afterRender = renderHook(() => useSpatialLoadTiming(deckRef, chart, describeLoad)).result.current;
    });

    /** Draw a frame of `deck` at time `at`. */
    function frame(at: number, deck: DeckGLRef) {
        now = at;
        deckRef.current = deck;
        afterRender();
    }

    it("reports the first load from when the chart opened, once", () => {
        frame(1800, deckAt(0, false));
        expect(recordSpatialLoaded).not.toHaveBeenCalled();

        frame(2500, deckAt(0, true));
        frame(2600, deckAt(0, true));

        expect(recordSpatialLoaded).toHaveBeenCalledTimes(1);
        expect(recordSpatialLoaded).toHaveBeenCalledWith("first", 1500, describeLoad());
    });

    it("does not report a second first load when the canvas remounts", () => {
        frame(2500, deckAt(0, true));
        recordSpatialLoaded.mockClear();

        afterRender = renderHook(() => useSpatialLoadTiming(deckRef, chart, describeLoad)).result.current;
        frame(10000, deckAt(0, false));
        frame(10500, deckAt(0, true));

        expect(recordSpatialLoaded).not.toHaveBeenCalled();
    });

    it("reports a settle from the last frame that moved, when tiles had to load", () => {
        frame(2000, deckAt(0, true));
        recordSpatialLoaded.mockClear();

        frame(3000, deckAt(1, true, false));
        frame(3100, deckAt(2, true, false));
        frame(3300, deckAt(2, true, false));
        frame(3400, deckAt(2, true, true));

        expect(recordSpatialLoaded).toHaveBeenCalledTimes(1);
        expect(recordSpatialLoaded).toHaveBeenCalledWith("settle", 300, describeLoad());
    });

    it("reports nothing for a move that needed no new tiles", () => {
        frame(2000, deckAt(0, true));
        recordSpatialLoaded.mockClear();

        frame(3000, deckAt(1, true));
        frame(3100, deckAt(1, true));

        expect(recordSpatialLoaded).not.toHaveBeenCalled();
    });
});
