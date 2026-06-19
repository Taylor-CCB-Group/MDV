import type { Layer } from "@deck.gl/core";
import type { SpatialData } from "@spatialdata/core";
import {
    RENDER_STACK_SCHEMA_VERSION,
    type RenderStack,
    type RenderStackEntry,
} from "@spatialdata/layers";
import { describe, expect, test, vi } from "vitest";

vi.mock("@spatialdata/vis", () => ({
    layerConfig: (
        type: string,
        base: Record<string, unknown>,
        props: Record<string, unknown> = {},
    ) => ({ ...props, ...base, type }),
    renderStackToLayerInputs: (stack: RenderStack) => ({
        layers: {},
        layerOrder: stack.entries
            .filter((entry) => entry.kind === "spatial")
            .map((entry) => entry.id),
    }),
}));

vi.mock("@/react/context", () => ({
    useChart: vi.fn(),
}));

vi.mock("@/react/hooks", () => ({
    useConfig: vi.fn(),
}));

import {
    createRenderStackLayerInputsCache,
    resolveCachedHostDeckLayers,
    syncRenderStackLayerInputs,
} from "@/react/spatialdata/render_stack_adapter";
import {
    createDefaultRenderStack,
    normalizeRenderStack,
} from "@/react/spatialdata/render_stack_defaults";
import {
    patchRenderStackEntry,
    removeRenderStackEntry,
} from "@/react/spatialdata/render_stack_control";
import { deckHostLayerId } from "@/react/spatialdata/host_overlay_ids";

type SpatialEntry = Extract<RenderStackEntry, { kind: "spatial" }>;
type HostEntry = Extract<RenderStackEntry, { kind: "host" }>;

function spatialEntry({
    id,
    elementKey,
    props = {},
    elementType = "image",
}: {
    id: string;
    elementKey: string;
    props?: Record<string, unknown>;
    elementType?: SpatialEntry["source"]["elementType"];
}): SpatialEntry {
    return {
        kind: "spatial",
        id,
        visible: true,
        source: { elementType, elementKey },
        props: { opacity: 1, ...props },
    };
}

function hostEntry(id: string, visible = true): HostEntry {
    return {
        kind: "host",
        id,
        visible,
        source: { hostLayerId: id },
        props: {},
    };
}

function renderStack(entries: RenderStackEntry[]): RenderStack {
    return {
        schemaVersion: RENDER_STACK_SCHEMA_VERSION,
        entries,
    };
}

function fakeDeckLayer(id: string) {
    const clonedLayer = { id: `${id}-clone` } as Layer;
    const clone = vi.fn(() => clonedLayer);
    const sourceLayer = { id, clone } as unknown as Layer;
    return { clonedLayer, clone, sourceLayer };
}

function fakeSpatialData(): SpatialData {
    return {
        images: {
            imageA: { coordinateSystems: ["global"] },
            imageB: { coordinateSystems: ["global"] },
        },
        shapes: {},
        points: {},
        labels: {},
    } as unknown as SpatialData;
}

describe("render stack adapter", () => {
    test("keeps layer object identity across cosmetic edits", () => {
        const imageEntry = spatialEntry({ id: "image-a", elementKey: "imageA" });
        const stack = renderStack([imageEntry]);
        const cache = createRenderStackLayerInputsCache();

        const first = syncRenderStackLayerInputs(stack, cache);
        const layers = first.layers;
        const layer = first.layers["image-a"];

        imageEntry.props.opacity = 0.25;
        const second = syncRenderStackLayerInputs(stack, cache);

        expect(second.layers).toBe(layers);
        expect(second.layers["image-a"]).toBe(layer);
        expect(second.layers["image-a"]?.opacity).toBe(0.25);
    });

    test("removes deleted spatial layers and updates layer order", () => {
        const imageA = spatialEntry({ id: "image-a", elementKey: "imageA" });
        const imageB = spatialEntry({ id: "image-b", elementKey: "imageB" });
        const stack = renderStack([imageA, imageB]);
        const cache = createRenderStackLayerInputsCache();

        expect(syncRenderStackLayerInputs(stack, cache).layerOrder).toEqual([
            "image-a",
            "image-b",
        ]);

        removeRenderStackEntry(stack, "image-a");
        const next = syncRenderStackLayerInputs(stack, cache);

        expect(next.layers["image-a"]).toBeUndefined();
        expect(next.layerOrder).toEqual(["image-b"]);
    });

    test("resolves only visible host entries and reuses host clones", () => {
        const visibleHostId = deckHostLayerId("scatter");
        const hiddenHostId = deckHostLayerId("selection");
        const stack = renderStack([
            hostEntry(visibleHostId),
            hostEntry(hiddenHostId, false),
        ]);
        const { clone, clonedLayer, sourceLayer } = fakeDeckLayer("scatter");
        const resolver = vi.fn(() => sourceLayer);

        const first = resolveCachedHostDeckLayers(stack, resolver);
        const second = resolveCachedHostDeckLayers(stack, resolver);

        expect(first).toEqual([clonedLayer]);
        expect(second).toEqual([clonedLayer]);
        expect(first[0]).toBe(second[0]);
        expect(clone).toHaveBeenCalledTimes(1);
        expect(resolver).toHaveBeenCalledTimes(2);
    });
});

describe("render stack control", () => {
    test("patches entry props without replacing unrelated entries", () => {
        const imageA = spatialEntry({ id: "image-a", elementKey: "imageA" });
        const imageB = spatialEntry({ id: "image-b", elementKey: "imageB" });
        const stack = renderStack([imageA, imageB]);

        patchRenderStackEntry(stack, "image-a", {
            visible: false,
            props: { opacity: 0.5 },
        });

        expect(stack.entries[0]).toBe(imageA);
        expect(stack.entries[1]).toBe(imageB);
        expect(imageA.visible).toBe(false);
        expect(imageA.props.opacity).toBe(0.5);
        expect(imageB.props.opacity).toBe(1);
    });
});

describe("render stack defaults", () => {
    test("normalizes defaults without replacing saved entries", () => {
        const savedImage = spatialEntry({
            id: "spatialdata-image-imageB",
            elementKey: "imageB",
            props: { opacity: 0.4 },
        });
        const savedStack = renderStack([savedImage]);

        const normalized = normalizeRenderStack(
            savedStack,
            fakeSpatialData(),
            "global",
        );
        const defaultStack = createDefaultRenderStack(fakeSpatialData(), "global");

        expect(normalized.entries[0]).toBe(savedImage);
        expect(normalized.entries.map((entry) => entry.id)).toEqual([
            "spatialdata-image-imageB",
            ...defaultStack.entries.map((entry) => entry.id),
        ]);
    });
});
