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
    renderStackOrder: (stack: RenderStack | undefined, fallback: string[]) =>
        stack
            ? stack.entries
                  .filter((entry) => entry.kind === "spatial")
                  .map((entry) => entry.id)
            : fallback,
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
    createHostOnlyRenderStack,
    normalizeRenderStack,
} from "@/react/spatialdata/render_stack_defaults";
import {
    patchRenderStackEntry,
    removeRenderStackEntry,
    seedRenderStackFromSpatialData,
} from "@/react/spatialdata/render_stack_control";
import { deckHostLayerId } from "@/react/spatialdata/host_overlay_ids";
import {
    renderStackSpatialRevision,
    touchRenderStackEntry,
} from "@/react/spatialdata/render_stack_observe";

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

    test("spatial revision changes when image props change in place", () => {
        const imageEntry = spatialEntry({
            id: "image-a",
            elementKey: "imageA",
            props: {
                channels: {
                    contrastLimits: [[0, 255]],
                    colors: [[255, 0, 0]],
                    channelsVisible: [true],
                    selections: [{ c: 0 }],
                },
            },
        });
        const stack = renderStack([imageEntry]);
        touchRenderStackEntry(imageEntry);
        const before = renderStackSpatialRevision(stack);

        patchRenderStackEntry(stack, "image-a", {
            props: {
                channels: {
                    contrastLimits: [[10, 200]],
                    colors: [[255, 0, 0]],
                    channelsVisible: [true],
                    selections: [{ c: 0 }],
                },
                vivLayerProps: { brightness: [0.5], contrast: [0.6] },
            },
        });
        touchRenderStackEntry(imageEntry);
        const afterChannels = renderStackSpatialRevision(stack);

        patchRenderStackEntry(stack, "image-a", {
            props: { vivLayerProps: { brightness: [0.7], contrast: [0.6] } },
        });
        touchRenderStackEntry(imageEntry);
        const afterTone = renderStackSpatialRevision(stack);

        expect(afterChannels).not.toBe(before);
        expect(afterTone).not.toBe(afterChannels);
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
    test("createHostOnlyRenderStack seeds host overlay entries", () => {
        const stack = createHostOnlyRenderStack();
        expect(stack.entries.length).toBeGreaterThan(0);
        expect(stack.entries.every((entry) => entry.kind === "host")).toBe(true);
    });

    test("seedRenderStackFromSpatialData adds the first image layer for brand-new charts", () => {
        const config = {
            type: "SpatialDataMdvRegionReact",
            renderStack: createHostOnlyRenderStack(),
        } as const;
        const chart = {
            seedDefaultSpatialLayers: true,
            bumpRenderStackGeneration: vi.fn(),
            finishDefaultSpatialLayerSeed: vi.fn(),
        };
        seedRenderStackFromSpatialData(
            config as never,
            fakeSpatialData(),
            "global",
            chart as never,
        );
        expect(
            config.renderStack?.entries.some(
                (entry) => entry.kind === "spatial" && entry.source.elementType === "image",
            ),
        ).toBe(true);
        expect(chart.finishDefaultSpatialLayerSeed).toHaveBeenCalled();
    });

    test("seedRenderStackFromSpatialData keeps persisted host-only stacks spatial-free", () => {
        const config = {
            type: "SpatialDataMdvRegionReact",
            renderStack: createHostOnlyRenderStack(),
        } as const;
        const chart = {
            seedDefaultSpatialLayers: false,
            bumpRenderStackGeneration: vi.fn(),
            finishDefaultSpatialLayerSeed: vi.fn(),
        };
        seedRenderStackFromSpatialData(
            config as never,
            fakeSpatialData(),
            "global",
            chart as never,
        );
        expect(
            config.renderStack?.entries.every((entry) => entry.kind === "host"),
        ).toBe(true);
        expect(chart.finishDefaultSpatialLayerSeed).not.toHaveBeenCalled();
    });

    test("seedRenderStackFromSpatialData keeps the same renderStack object", () => {
        const config = {
            type: "SpatialDataMdvRegionReact",
            renderStack: createHostOnlyRenderStack(),
        } as const;
        const stackBefore = config.renderStack;
        const chart = {
            seedDefaultSpatialLayers: true,
            bumpRenderStackGeneration: vi.fn(),
            finishDefaultSpatialLayerSeed: vi.fn(),
        };
        seedRenderStackFromSpatialData(
            config as never,
            fakeSpatialData(),
            "global",
            chart as never,
        );
        expect(config.renderStack).toBe(stackBefore);
    });

    test("normalizeRenderStack adds missing host overlays without injecting spatial layers", () => {
        const savedImage = spatialEntry({
            id: "spatialdata-image-imageB",
            elementKey: "imageB",
            props: { opacity: 0.4 },
        });
        const savedStack = renderStack([savedImage]);

        const normalized = normalizeRenderStack(savedStack);

        expect(normalized.entries[0]).toBe(savedImage);
        expect(normalized.entries.some((entry) => entry.kind === "host")).toBe(true);
        expect(
            normalized.entries.filter((entry) => entry.kind === "spatial").length,
        ).toBe(1);
    });
});
