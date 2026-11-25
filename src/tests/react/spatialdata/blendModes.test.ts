import type { Layer } from "@deck.gl/core";
import { ScatterplotLayer } from "@deck.gl/layers";
import { IsolatedGroupLayer, type RenderStack, type RenderStackEntry } from "@spatialdata/layers";
import { renderStackToDeckLayers } from "@spatialdata/vis";
import { observable } from "mobx";
import { describe, expect, test } from "vitest";

import {
    compositingBadge,
    isInsideIsolatedGroup,
    leafBlendMode,
    renderStackAncestors,
    setLeafBlendMode,
} from "@/react/spatialdata/blend_modes";
import { DATA_POINTS_GROUP_ID, deckHostLayerId } from "@/react/spatialdata/host_overlay_ids";
import { applyIsolatedGroupFix, BLEND_WRAPPER_SUFFIX } from "@/react/spatialdata/isolated_group_fix";
import { renderStackStructure } from "@/react/spatialdata/render_stack_adapter";

function spatial(id: string, props: Record<string, unknown> = {}): RenderStackEntry {
    return { kind: "spatial", id, visible: true, source: { elementType: "shapes", elementKey: id }, props };
}

function group(id: string, props: Record<string, unknown>, children: RenderStackEntry[]): RenderStackEntry {
    return { kind: "group", id, visible: true, props, children };
}

function host(id: string): RenderStackEntry {
    return { kind: "host", id, visible: true, source: { hostLayerId: id }, props: {} };
}

/** image, then an isolated group holding a pass-through group holding `mask`, then `loose` at the top. */
function stack(maskProps: Record<string, unknown> = {}): RenderStack {
    return {
        schemaVersion: 2,
        entries: [
            spatial("image"),
            group("outer", { isolate: true }, [spatial("cells"), group("inner", {}, [spatial("mask", maskProps)])]),
            spatial("loose", { blendMode: "multiply" }),
        ],
    };
}

describe("isInsideIsolatedGroup", () => {
    test("is true under an isolated ancestor at any depth, through pass-through groups", () => {
        const s = stack();
        expect(renderStackAncestors(s, "mask")?.map((g) => g.id)).toEqual(["outer", "inner"]);
        expect(isInsideIsolatedGroup(s, "mask")).toBe(true);
        expect(isInsideIsolatedGroup(s, "inner")).toBe(true);
        expect(isInsideIsolatedGroup(s, "outer")).toBe(false);
        expect(isInsideIsolatedGroup(s, "loose")).toBe(false);
        expect(renderStackAncestors(s, "missing")).toBeUndefined();
    });

    test("a blend mode other than pass-through isolates the group", () => {
        const s: RenderStack = {
            schemaVersion: 2,
            entries: [group("g", { blendMode: "screen" }, [spatial("a")])],
        };
        expect(isInsideIsolatedGroup(s, "a")).toBe(true);
    });
});

describe("setLeafBlendMode", () => {
    test("stores normal as no mode, on plain and observable stacks", () => {
        for (const s of [stack(), observable(stack())]) {
            expect(setLeafBlendMode(s, "mask", "cut-out")).toBe(true);
            const outer = s.entries[1];
            const inner = outer.kind === "group" ? outer.children[1] : undefined;
            const mask = inner?.kind === "group" ? inner.children[0] : undefined;
            expect(mask && leafBlendMode(mask)).toBe("cut-out");
            expect(setLeafBlendMode(s, "mask", "cut-out")).toBe(false);
            expect(setLeafBlendMode(s, "mask", "normal")).toBe(true);
            expect(mask && "blendMode" in mask.props).toBe(false);
            expect(setLeafBlendMode(s, "outer", "multiply")).toBe(false);
        }
    });
});

describe("compositingBadge", () => {
    test("shows group modes, isolation and leaf modes, marking those drawn as normal", () => {
        const s = stack({ blendMode: "clip" });
        const [, outer, loose] = s.entries;
        expect(compositingBadge(s, outer)).toMatchObject({ label: "One layer", inactive: false });
        expect(compositingBadge(s, loose)).toMatchObject({ label: "Multiply", inactive: true });
        const mask =
            outer.kind === "group" && outer.children[1].kind === "group" ? outer.children[1].children[0] : loose;
        expect(compositingBadge(s, mask)).toMatchObject({ label: "Clip", inactive: false });
        expect(compositingBadge(s, s.entries[0])).toBeNull();
    });

    test("a top-level group's operator composites as normal", () => {
        const s: RenderStack = { schemaVersion: 2, entries: [group("g", { blendMode: "intersect" }, [spatial("a")])] };
        expect(compositingBadge(s, s.entries[0])).toMatchObject({ label: "Intersect", inactive: true });
    });

    test("the data points row shows its group's compositing", () => {
        const scatter = host(deckHostLayerId("scatter"));
        const s: RenderStack = {
            schemaVersion: 2,
            entries: [group(DATA_POINTS_GROUP_ID, { blendMode: "multiply", opacity: 0.5 }, [scatter])],
        };
        expect(compositingBadge(s, scatter)).toMatchObject({ label: "Multiply", inactive: false });
    });
});

/**
 * Which deck layer class a layer is, by its `layerName`: under vitest's dependency optimizer,
 * @spatialdata/vis can hold a copy of @spatialdata/layers of its own, so `instanceof` against
 * the class imported here can fail for a layer vis built.
 */
function layerName(layer: unknown): unknown {
    return layer instanceof Object ? Reflect.get(layer.constructor, "layerName") : undefined;
}

describe("renderStackStructure", () => {
    test("keeps a layer's blend mode and drops its other props", () => {
        const s = stack({ blendMode: "cut-out", fillColor: [1, 2, 3], opacity: 0.5 });
        const structure = renderStackStructure(s);
        const outer = structure.entries[1];
        const inner = outer.kind === "group" ? outer.children[1] : undefined;
        const mask = inner?.kind === "group" ? inner.children[0] : undefined;
        expect(mask?.props).toEqual({ blendMode: "cut-out" });
        expect(structure.entries[0].props).toEqual({});
    });

    test("sd.js draws a blended layer inside an isolated group as a one-member group with its mode", () => {
        const structure = renderStackStructure(stack({ blendMode: "cut-out" }));
        const entryLayers = new Map<string, Layer[]>(
            ["image", "cells", "mask", "loose"].map((id) => [id, [new ScatterplotLayer({ id, data: [] })]]),
        );
        const layers = renderStackToDeckLayers(structure, { entryLayers });
        expect(layers.map((layer) => layer.id)).toEqual(["image", "outer", "loose"]);
        // outside any isolated group, the multiply on `loose` is not drawn
        expect(layers[2]).toBeInstanceOf(ScatterplotLayer);
        const outer = layers[1];
        expect(layerName(outer)).toBe(IsolatedGroupLayer.layerName);
        const members = Reflect.get(outer.props, "layers");
        const mask: unknown = Array.isArray(members) ? members[1] : undefined;
        expect(layerName(mask)).toBe(IsolatedGroupLayer.layerName);
        expect(mask instanceof Object && Reflect.get(Reflect.get(mask, "props"), "blendMode")).toBe("cut-out");
    });
});

describe("applyIsolatedGroupFix", () => {
    test("gives a blended layer's one-member wrapper a sublayer id apart from the layer's", () => {
        applyIsolatedGroupFix();
        const member = new ScatterplotLayer({ id: "shapes", data: [] });
        const wrapper = new IsolatedGroupLayer({ id: "shapes", layers: [member], blendMode: "cut-out" });
        const plain = new ScatterplotLayer({ id: "image", data: [] });
        const outer = new IsolatedGroupLayer({ id: "outer", layers: [plain, wrapper] });
        const sublayers = outer.renderLayers();
        const ids: unknown[] = [];
        const collect = (item: unknown): void => {
            if (Array.isArray(item)) item.forEach(collect);
            else if (item) ids.push(Reflect.get(Object(item), "id"));
        };
        collect(sublayers);
        expect(ids).toEqual(["outer/image", `outer/shapes${BLEND_WRAPPER_SUFFIX}`, "outer-composite"]);
    });

    test("never gives the wrapper a sibling's id", () => {
        applyIsolatedGroupFix();
        const wrapper = new IsolatedGroupLayer({
            id: "shapes",
            layers: [new ScatterplotLayer({ id: "shapes", data: [] })],
            blendMode: "multiply",
        });
        const sibling = new ScatterplotLayer({ id: "shapes-blend", data: [] });
        const outer = new IsolatedGroupLayer({ id: "outer", layers: [sibling, wrapper] });
        const ids: unknown[] = [];
        const collect = (item: unknown): void => {
            if (Array.isArray(item)) item.forEach(collect);
            else if (item) ids.push(Reflect.get(Object(item), "id"));
        };
        collect(outer.renderLayers());
        expect(ids).toEqual(["outer/shapes-blend", "outer/shapes-blend-2", "outer-composite"]);
        expect(new Set(ids).size).toBe(ids.length);
    });
});
