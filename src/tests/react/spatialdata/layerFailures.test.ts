import {
    CAUSE_WINDOW_MS,
    type LayerFailureSink,
    LayerFailureTracker,
    entryForLabel,
} from "@/react/spatialdata/layer_failures";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

// Messages as luma reports them for a points layer whose palette texture is too wide.
const POINTS = "spatialdata-points-transcripts-#L1Htwz";
const TOO_BIG =
    "Texture size ([Extent3D width:12448, height:1, depthOrArrayLayers:1]) exceeded maximum texture size " +
    '([Extent3D width:8192, height:8192, depthOrArrayLayers:256]).\n - While validating [TextureDescriptor ""Texture-11""].';
const CONSEQUENCE = '[Invalid Texture "Texture-11"] is invalid due to a previous error.\n - While calling CreateView.';
const BIND_GROUP =
    '[Invalid TextureView "Texture-12"] is invalid due to a previous error.\n - While validating [BindGroupDescriptor ""' +
    `${POINTS}-#spatial-detail-b3y4whn3l#-scatter-cached/group0[scatterplot,project]""]`;
const DRAW = `Debug group stack:\n > "Model(${POINTS}-#spatial-detail-b3y4whn3l#-scatter).draw(RenderPass:"RenderPass-17") "`;

const ENTRY_IDS = ["spatialdata-image-morphology_focus-#LF3nxE", POINTS, "deck:grey_scatter", "deck:scatter"];

describe("entryForLabel", () => {
    it("finds the entry in luma resource labels", () => {
        expect(entryForLabel(BIND_GROUP, ENTRY_IDS)).toBe(POINTS);
        expect(entryForLabel(DRAW, ENTRY_IDS)).toBe(POINTS);
    });

    it("finds the entry from a deck layer id, including sublayers", () => {
        expect(entryForLabel("deck:scatter-#spatial-detail-x#-spatial.scatterplot", ENTRY_IDS)).toBe("deck:scatter");
        expect(entryForLabel("deck:grey_scatter-#spatial-detail-x#", ENTRY_IDS)).toBe("deck:grey_scatter");
    });

    it("needs a whole label prefix, not a substring", () => {
        expect(entryForLabel("mydeck:scatterplot", ENTRY_IDS)).toBeUndefined();
        expect(entryForLabel(TOO_BIG, ENTRY_IDS)).toBeUndefined();
    });
});

describe("LayerFailureTracker", () => {
    let fail: ReturnType<typeof vi.fn<LayerFailureSink["fail"]>>;
    let canvas: ReturnType<typeof vi.fn<LayerFailureSink["canvas"]>>;
    let log: ReturnType<typeof vi.fn<(error: Error) => void>>;
    let now = 0;
    let tracker: LayerFailureTracker;

    beforeEach(() => {
        vi.useFakeTimers();
        now = 0;
        fail = vi.fn<LayerFailureSink["fail"]>();
        canvas = vi.fn<LayerFailureSink["canvas"]>();
        log = vi.fn<(error: Error) => void>();
        tracker = new LayerFailureTracker({ entryIds: () => ENTRY_IDS, fail, canvas }, () => now, log);
    });
    afterEach(() => {
        tracker.dispose();
        vi.useRealTimers();
    });

    it("blames the layer for the unattributed error that came first", () => {
        tracker.device(new Error(TOO_BIG));
        tracker.device(new Error(CONSEQUENCE));
        now = 5;
        tracker.device(new Error(BIND_GROUP));
        expect(fail).toHaveBeenCalledWith(POINTS, { reason: TOO_BIG, detail: BIND_GROUP, source: "webgpu" });
        vi.advanceTimersByTime(CAUSE_WINDOW_MS * 2);
        expect(canvas).not.toHaveBeenCalled();
    });

    it("reports an error nothing claims on the canvas", () => {
        tracker.device(new Error(TOO_BIG));
        vi.advanceTimersByTime(CAUSE_WINDOW_MS);
        expect(canvas).toHaveBeenCalledWith(TOO_BIG);
        expect(fail).not.toHaveBeenCalled();
    });

    it("does not report consequences on their own", () => {
        tracker.device(new Error(CONSEQUENCE));
        vi.advanceTimersByTime(CAUSE_WINDOW_MS * 2);
        expect(canvas).not.toHaveBeenCalled();
    });

    it("uses the attributed message alone when no cause is held", () => {
        tracker.device(new Error(DRAW));
        expect(fail).toHaveBeenCalledWith(POINTS, { reason: DRAW, source: "webgpu" });
    });

    it("attributes deck layer errors by layer id, or to the canvas, keeping their stacks", () => {
        const shaderError = new Error("shader failed");
        tracker.layer(shaderError, `${POINTS}-#spatial-detail-x#-scatter`);
        expect(fail).toHaveBeenCalledWith(POINTS, { reason: "shader failed", source: "deck", stack: shaderError.stack });
        const unattributed = new Error("no layer");
        tracker.layer(unattributed, undefined);
        expect(canvas).toHaveBeenCalledWith("no layer", unattributed.stack);
    });

    it("logs each kind of message once, however often it repeats", () => {
        for (let frame = 0; frame < 100; frame++) tracker.device(new Error(DRAW.replace("17", String(frame))));
        expect(log).toHaveBeenCalledTimes(1);
    });
});
