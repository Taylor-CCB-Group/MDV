import type { Layer } from "@deck.gl/core";
import type { SpatialData } from "@spatialdata/core";
import {
    getRenderStackEntryIds,
    RENDER_STACK_SCHEMA_VERSION,
    type RenderStack,
    type RenderStackEntry,
    type RenderStackV1,
} from "@spatialdata/layers";
import { describe, expect, test, vi } from "vitest";

vi.mock("@/react/context", () => ({
    useChart: vi.fn(),
}));

vi.mock("@/react/hooks", () => ({
    useConfig: vi.fn(),
}));

import { DECK_OVERLAY_IDS, deckHostLayerId } from "@/react/spatialdata/host_overlay_ids";
import { flattenLayerTree, projectLayerTreeDrop } from "@/react/spatialdata/layer_tree";
import {
    createRenderStackLayerInputsCache,
    resolveCachedHostDeckLayers,
    syncRenderStackLayerInputs,
} from "@/react/spatialdata/render_stack_adapter";
import {
    DATA_POINTS_GROUP_ID,
    DEFAULT_COMPOSITING,
    insertRenderStackGroup,
    insertSpatialRenderStackEntry,
    moveRenderStackEntry,
    patchRenderStackEntry,
    patchRenderStackEntryPropsInPlace,
    removeRenderStackEntry,
    seedRenderStackFromSpatialData,
    setDataPointsCompositing,
    ungroupRenderStackEntry,
} from "@/react/spatialdata/render_stack_control";
import {
    createHostOnlyRenderStack,
    mergeHostOverlayEntriesInPlace,
    migrateSavedRenderStack,
    normalizeRenderStack,
    syncOverlayFollowerVisibility,
} from "@/react/spatialdata/render_stack_defaults";
import {
    customDisplayName,
    numberedDefaultName,
    renderStackEntryDisplayName,
} from "@/react/spatialdata/render_stack_display";
import { renderStackSpatialRevision, touchRenderStackEntry } from "@/react/spatialdata/render_stack_observe";

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
    let cloneCount = 0;
    const clone = vi.fn((props?: { id?: string }) => {
        cloneCount += 1;
        return { id: props?.id ?? `${id}-clone-${cloneCount}` } as Layer;
    });
    const sourceLayer = { id, clone } as unknown as Layer;
    return { clone, sourceLayer };
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

    test("replaces a layer config when table-driven annotation props change", () => {
        const labelsEntry = spatialEntry({
            id: "labels-a",
            elementKey: "cell_labels",
            elementType: "labels",
            props: {
                fillColorByColumn: {
                    columnName: "Expressed_genes",
                    mode: "categorical",
                },
            },
        });
        const stack = renderStack([labelsEntry]);
        const cache = createRenderStackLayerInputsCache();

        const first = syncRenderStackLayerInputs(stack, cache);
        const layers = first.layers;
        const layer = first.layers["labels-a"];

        patchRenderStackEntry(stack, "labels-a", {
            props: {
                fillColorByColumn: {
                    columnName: "Leiden",
                    mode: "categorical",
                },
            },
        });
        const second = syncRenderStackLayerInputs(stack, cache);

        expect(second.layers).toBe(layers);
        expect(second.layers["labels-a"]).not.toBe(layer);
        expect(second.layers["labels-a"]).toMatchObject({
            fillColorByColumn: {
                columnName: "Leiden",
                mode: "categorical",
            },
        });
    });

    test("removes deleted spatial layers and updates layer order", () => {
        const imageA = spatialEntry({ id: "image-a", elementKey: "imageA" });
        const imageB = spatialEntry({ id: "image-b", elementKey: "imageB" });
        const stack = renderStack([imageA, imageB]);
        const cache = createRenderStackLayerInputsCache();

        expect(syncRenderStackLayerInputs(stack, cache).layerOrder).toEqual(["image-a", "image-b"]);

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

    test("spatial revision changes when non-image layer props change in place", () => {
        const labelsEntry = spatialEntry({
            id: "labels-a",
            elementKey: "cell_labels",
            elementType: "labels",
        });
        const stack = renderStack([labelsEntry]);
        touchRenderStackEntry(labelsEntry);
        const before = renderStackSpatialRevision(stack);

        patchRenderStackEntry(stack, "labels-a", {
            props: {
                fillColorByColumn: {
                    columnName: "Expressed_genes",
                    mode: "categorical",
                },
            },
        });
        touchRenderStackEntry(labelsEntry);
        const afterFillColumn = renderStackSpatialRevision(stack);

        patchRenderStackEntry(stack, "labels-a", {
            props: { tooltipFields: ["Expressed_genes", "Leiden"] },
        });
        touchRenderStackEntry(labelsEntry);
        const afterTooltipFields = renderStackSpatialRevision(stack);

        expect(afterFillColumn).not.toBe(before);
        expect(afterTooltipFields).not.toBe(afterFillColumn);
    });

    test("patching a layer prop to undefined clears it and changes spatial revision", () => {
        const labelsEntry = spatialEntry({
            id: "labels-a",
            elementKey: "cell_labels",
            elementType: "labels",
            props: {
                fillColorByColumn: {
                    columnName: "Expressed_genes",
                    mode: "categorical",
                },
            },
        });
        const stack = renderStack([labelsEntry]);
        touchRenderStackEntry(labelsEntry);
        const before = renderStackSpatialRevision(stack);

        patchRenderStackEntry(stack, "labels-a", {
            props: { fillColorByColumn: undefined },
        });
        touchRenderStackEntry(labelsEntry);
        const after = renderStackSpatialRevision(stack);

        expect(labelsEntry.props.fillColorByColumn).toBeUndefined();
        expect(after).not.toBe(before);
    });

    test("resolves only visible host entries and refreshes host clones", () => {
        const visibleHostId = deckHostLayerId("scatter");
        const hiddenHostId = deckHostLayerId("selection");
        const stack = renderStack([hostEntry(visibleHostId), hostEntry(hiddenHostId, false)]);
        const { clone, sourceLayer } = fakeDeckLayer("scatter");
        const resolver = vi.fn(() => sourceLayer);

        const first = resolveCachedHostDeckLayers(stack, resolver);
        const second = resolveCachedHostDeckLayers(stack, resolver);

        expect(first).toHaveLength(1);
        expect(second).toHaveLength(1);
        expect(first[0]).not.toBe(second[0]);
        expect(first[0]?.id).toBe(visibleHostId);
        expect(second[0]?.id).toBe(visibleHostId);
        expect(clone).toHaveBeenCalledTimes(2);
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

    test("merges labels channel appearance into saved channels", () => {
        const props: Record<string, unknown> = { channels: { colors: [[255, 0, 0]], channelOpacities: [0.2] } };
        patchRenderStackEntryPropsInPlace(props, {
            channels: { channelsFilled: [false], channelStrokeWidths: [3] },
        });
        expect(props.channels).toEqual({
            colors: [[255, 0, 0]],
            channelOpacities: [0.2],
            channelsFilled: [false],
            channelStrokeWidths: [3],
        });
    });
});

describe("render stack display names", () => {
    test("a user-chosen name wins; blank names fall back to the element", () => {
        const named = spatialEntry({
            id: "a",
            elementKey: "cells",
            elementType: "shapes",
            props: { displayName: " Nuclei " },
        });
        const blank = spatialEntry({
            id: "b",
            elementKey: "cells",
            elementType: "shapes",
            props: { displayName: "  " },
        });
        expect(renderStackEntryDisplayName(named)).toBe("Nuclei");
        expect(customDisplayName(blank)).toBeUndefined();
        expect(renderStackEntryDisplayName(blank)).toBe("shapes: cells");
    });

    test("unnamed entries of the same element are numbered in stack order", () => {
        const first = spatialEntry({ id: "a", elementKey: "cells", elementType: "shapes" });
        const second = spatialEntry({ id: "b", elementKey: "cells", elementType: "shapes" });
        const other = spatialEntry({ id: "c", elementKey: "imageA" });
        const stack = renderStack([first, other, second]);
        expect(renderStackEntryDisplayName(first, stack)).toBe("shapes: cells (1)");
        expect(renderStackEntryDisplayName(second, stack)).toBe("shapes: cells (2)");
        expect(renderStackEntryDisplayName(other, stack)).toBe("image: imageA");

        second.props.displayName = "Outlines";
        expect(renderStackEntryDisplayName(second, stack)).toBe("Outlines");
        // the default stays numbered, e.g. as the placeholder of a named entry
        expect(numberedDefaultName(second, stack)).toBe("shapes: cells (2)");
    });

    test("the name does not reach the viewer's layer config", () => {
        const entry = spatialEntry({ id: "a", elementKey: "imageA", props: { displayName: "DAPI" } });
        const inputs = syncRenderStackLayerInputs(renderStack([entry]), createRenderStackLayerInputsCache());
        expect(inputs.layers.a).toBeDefined();
        expect("displayName" in (inputs.layers.a ?? {})).toBe(false);
    });
});

describe("render stack defaults", () => {
    test("a restored stack's grey scatter follows the data points' visibility", () => {
        const scatter = hostEntry(deckHostLayerId("scatter"), false);
        const grey = hostEntry(deckHostLayerId("grey_scatter"), true);
        const gates = hostEntry(deckHostLayerId("gate_display"), false);
        syncOverlayFollowerVisibility(renderStack([grey, scatter, gates]));
        expect(grey.visible).toBe(false);
        expect(gates.visible).toBe(false);
        scatter.visible = true;
        syncOverlayFollowerVisibility(renderStack([grey, scatter, gates]));
        expect(grey.visible).toBe(true);
        expect(gates.visible).toBe(false);
    });

    test("a grey scatter added to a saved stack follows the data points' visibility", () => {
        const scatter = hostEntry(deckHostLayerId("scatter"), false);
        const stack = renderStack([scatter]);
        expect(mergeHostOverlayEntriesInPlace(stack)).toBe(true);
        const grey = stack.entries.find((entry) => entry.id === deckHostLayerId("grey_scatter"));
        expect(grey?.visible).toBe(false);
        // Other added overlays keep their default.
        const gates = stack.entries.find((entry) => entry.id === deckHostLayerId("gate_display"));
        expect(gates?.visible).toBe(true);
    });

    test("an existing grey scatter keeps its saved visibility when merging", () => {
        const scatter = hostEntry(deckHostLayerId("scatter"), false);
        const grey = hostEntry(deckHostLayerId("grey_scatter"), true);
        const stack = renderStack([grey, scatter]);
        mergeHostOverlayEntriesInPlace(stack);
        expect(grey.visible).toBe(true);
    });

    test("a saved stack without density gets it just above the data points", () => {
        // A stack saved before density was its own overlay: an image under the MDV overlays.
        const image = {
            kind: "spatial",
            id: "image-1",
            visible: true,
            source: { elementType: "image", elementKey: "img" },
            props: {},
        } satisfies RenderStackEntry;
        const saved = ["grey_scatter", "scatter", "gate_display", "selection", "gate_labels"] as const;
        const stack = renderStack([image, ...saved.map((deckId) => hostEntry(deckHostLayerId(deckId)))]);
        expect(mergeHostOverlayEntriesInPlace(stack)).toBe(true);
        expect(stack.entries.map((entry) => entry.id)).toEqual([
            "image-1",
            "deck:grey_scatter",
            "deck:scatter",
            "deck:density",
            "deck:gate_display",
            "deck:selection",
            "deck:gate_labels",
        ]);
    });

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
            bumpRenderStackPropsGeneration: vi.fn(),
            finishDefaultSpatialLayerSeed: vi.fn(),
        };
        seedRenderStackFromSpatialData(config as never, fakeSpatialData(), "global", chart as never);
        expect(
            config.renderStack?.entries.some(
                (entry) => entry.kind === "spatial" && entry.source.elementType === "image",
            ),
        ).toBe(true);
        expect(chart.finishDefaultSpatialLayerSeed).toHaveBeenCalled();
    });

    test("seedRenderStackFromSpatialData gives the first image layer the chart's default image props", () => {
        const config = {
            type: "SpatialDataMdvRegionReact",
            renderStack: createHostOnlyRenderStack(),
        } as const;
        const defaultImageLayerProps = {
            channels: { colors: [[255, 0, 0]], contrastLimits: [[0, 100]] },
            vivLayerProps: { brightness: [0.6] },
        };
        const chart = {
            seedDefaultSpatialLayers: true,
            defaultImageLayerProps,
            bumpRenderStackGeneration: vi.fn(),
            bumpRenderStackPropsGeneration: vi.fn(),
            finishDefaultSpatialLayerSeed: vi.fn(),
        };
        seedRenderStackFromSpatialData(config as never, fakeSpatialData(), "global", chart as never);
        const image = config.renderStack.entries.find(
            (entry) => entry.kind === "spatial" && entry.source.elementType === "image",
        );
        expect(image?.props).toEqual({ opacity: 1, ...defaultImageLayerProps });
        // a copy, so editing the layer leaves the chart's props alone
        expect(image?.props.channels).not.toBe(defaultImageLayerProps.channels);
    });

    test("seedRenderStackFromSpatialData keeps persisted host-only stacks spatial-free", () => {
        const config = {
            type: "SpatialDataMdvRegionReact",
            renderStack: createHostOnlyRenderStack(),
        } as const;
        const chart = {
            seedDefaultSpatialLayers: false,
            bumpRenderStackGeneration: vi.fn(),
            bumpRenderStackPropsGeneration: vi.fn(),
            finishDefaultSpatialLayerSeed: vi.fn(),
        };
        seedRenderStackFromSpatialData(config as never, fakeSpatialData(), "global", chart as never);
        expect(config.renderStack?.entries.every((entry) => entry.kind === "host")).toBe(true);
        expect(chart.finishDefaultSpatialLayerSeed).not.toHaveBeenCalled();
        expect(chart.bumpRenderStackPropsGeneration).toHaveBeenCalledTimes(1);
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
            bumpRenderStackPropsGeneration: vi.fn(),
            finishDefaultSpatialLayerSeed: vi.fn(),
        };
        seedRenderStackFromSpatialData(config as never, fakeSpatialData(), "global", chart as never);
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
        expect(normalized.entries.filter((entry) => entry.kind === "spatial").length).toBe(1);
    });
});

describe("render stack groups", () => {
    const hostIds = DECK_OVERLAY_IDS.map(deckHostLayerId);

    function groupedStack(groupProps: Record<string, unknown> = {}, visible = true): RenderStack {
        return renderStack([
            spatialEntry({ id: "image-a", elementKey: "imageA" }),
            {
                kind: "group",
                id: "points",
                visible,
                props: groupProps,
                children: [
                    hostEntry(deckHostLayerId("grey_scatter")),
                    hostEntry(deckHostLayerId("scatter")),
                    spatialEntry({ id: "shapes-a", elementKey: "cells", elementType: "shapes" }),
                ],
            },
            ...hostIds.slice(2).map((id) => hostEntry(id)),
        ]);
    }

    test("a saved v1 stack loads unchanged, at v2", () => {
        // MDV's default stack as saved before v2: an image, then the host overlays.
        const v1: RenderStackV1 = {
            schemaVersion: 1,
            entries: [
                spatialEntry({ id: "image-a", elementKey: "imageA", props: { vivLayerProps: { contrast: 2 } } }),
                ...hostIds.map((id) => hostEntry(id, id !== deckHostLayerId("grey_scatter"))),
            ],
        };
        const before = structuredClone(v1);

        const loaded = migrateSavedRenderStack(v1);

        expect(loaded).toEqual({ ...before, schemaVersion: 2 });
        expect(normalizeRenderStack(loaded)).toEqual(loaded);
    });

    test("normalizing finds host entries inside groups, and migrates a v1 group", () => {
        const stack = groupedStack();
        expect(getRenderStackEntryIds(normalizeRenderStack(stack))).toEqual(getRenderStackEntryIds(stack));

        const v1 = {
            schemaVersion: 1,
            entries: [
                hostEntry(deckHostLayerId("grey_scatter")),
                { kind: "group", id: "g", visible: true, props: {}, children: [deckHostLayerId("grey_scatter")] },
            ],
        } satisfies RenderStackV1;
        const live = v1 as unknown as RenderStack;
        expect(mergeHostOverlayEntriesInPlace(live)).toBe(true);
        expect(live.schemaVersion).toBe(2);
        expect(getRenderStackEntryIds(live)).toEqual(["g", ...hostIds]);
    });

    test("a missing overlay goes just above the group holding the overlay it follows", () => {
        const stack = groupedStack();
        const group = stack.entries[1];
        stack.entries.splice(2, 1); // drop deck:density, which sits right after the group
        expect(stack.entries.map((entry) => entry.id)).not.toContain(deckHostLayerId("density"));

        expect(mergeHostOverlayEntriesInPlace(stack)).toBe(true);
        expect(stack.entries.map((entry) => entry.id)).toEqual([
            "image-a",
            "points",
            deckHostLayerId("density"),
            ...hostIds.slice(3),
        ]);
        expect(stack.entries[1]).toBe(group);
    });

    test("syncs spatial entries and resolves host entries at any depth, inheriting visibility", () => {
        const visible = syncRenderStackLayerInputs(groupedStack(), createRenderStackLayerInputsCache());
        expect(Object.keys(visible.layers)).toEqual(["image-a", "shapes-a"]);
        expect(visible.layerOrder).toEqual(["image-a", ...hostIds.slice(0, 2), "shapes-a", ...hostIds.slice(2)]);

        const hidden = groupedStack({}, false);
        const inputs = syncRenderStackLayerInputs(hidden, createRenderStackLayerInputsCache());
        expect(inputs.layers["shapes-a"]?.visible).toBe(false);

        const resolver = (entry: HostEntry) => fakeDeckLayer(entry.id).sourceLayer;
        expect(resolveCachedHostDeckLayers(groupedStack(), resolver).map((layer) => layer.id)).toEqual(hostIds);
        expect(resolveCachedHostDeckLayers(hidden, resolver).map((layer) => layer.id)).toEqual(hostIds.slice(2));
    });

    test("the data points' group exists only while their compositing is not the default", () => {
        const stack = createHostOnlyRenderStack();
        const isolated = { ...DEFAULT_COMPOSITING, opacity: 0.5, isolate: true };

        expect(setDataPointsCompositing(stack, isolated)).toBe(true);
        expect(getRenderStackEntryIds(stack)).toEqual([DATA_POINTS_GROUP_ID, ...hostIds]);
        expect(setDataPointsCompositing(stack, { ...isolated, blendMode: "multiply" })).toBe(false);
        expect(stack.entries[0]).toMatchObject({ props: { blendMode: "multiply" } });

        expect(setDataPointsCompositing(stack, DEFAULT_COMPOSITING)).toBe(true);
        expect(getRenderStackEntryIds(stack)).toEqual(hostIds);
        expect(setDataPointsCompositing(stack, DEFAULT_COMPOSITING)).toBe(false);
    });

    test("groups are inserted below the overlays, filled by moving entries in, and ungrouped in place", () => {
        const stack = renderStack([
            spatialEntry({ id: "image-a", elementKey: "imageA" }),
            spatialEntry({ id: "shapes-a", elementKey: "cells", elementType: "shapes" }),
            hostEntry(deckHostLayerId("scatter")),
        ]);
        const groupId = insertRenderStackGroup(stack);
        expect(stack.entries.map((entry) => entry.id)).toEqual(["image-a", "shapes-a", groupId, "deck:scatter"]);

        expect(moveRenderStackEntry(stack, "shapes-a", groupId, null)).toBe(true);
        expect(moveRenderStackEntry(stack, "image-a", groupId, "shapes-a")).toBe(true);
        expect(getRenderStackEntryIds(stack)).toEqual([groupId, "shapes-a", "image-a", "deck:scatter"]);
        expect(moveRenderStackEntry(stack, groupId, groupId, null)).toBe(false);

        removeRenderStackEntry(stack, "image-a");
        expect(getRenderStackEntryIds(stack)).toEqual([groupId, "shapes-a", "deck:scatter"]);

        expect(ungroupRenderStackEntry(stack, groupId)).toBe(true);
        expect(stack.entries.map((entry) => entry.id)).toEqual(["shapes-a", "deck:scatter"]);
    });

    test("new layers go on top of a group, or below the overlays and the data points' group", () => {
        const stack = createHostOnlyRenderStack();
        setDataPointsCompositing(stack, { ...DEFAULT_COMPOSITING, isolate: true });
        const groupId = insertRenderStackGroup(stack);
        expect(stack.entries.map((entry) => entry.id).slice(0, 2)).toEqual([groupId, DATA_POINTS_GROUP_ID]);

        const imageId = insertSpatialRenderStackEntry(stack, "image", "imageA", {}, groupId);
        const shapesId = insertSpatialRenderStackEntry(stack, "shapes", "cells", {}, groupId);
        expect(getRenderStackEntryIds(stack).slice(0, 3)).toEqual([groupId, imageId, shapesId]);
        const labelsId = insertSpatialRenderStackEntry(stack, "labels", "mask");
        expect(stack.entries.map((entry) => entry.id).slice(0, 3)).toEqual([groupId, labelsId, DATA_POINTS_GROUP_ID]);
    });

    test("a dragged row lands at the depth the rows around it allow", () => {
        const stack = renderStack([
            spatialEntry({ id: "a", elementKey: "a" }),
            {
                kind: "group",
                id: "g",
                visible: true,
                props: {},
                children: [spatialEntry({ id: "b", elementKey: "b" })],
            },
            spatialEntry({ id: "c", elementKey: "c" }),
        ]);
        const rows = flattenLayerTree(stack.entries, new Set());
        expect(rows.map((row) => `${row.id}@${row.depth}`)).toEqual(["a@0", "g@0", "b@1", "c@0"]);
        expect(flattenLayerTree(stack.entries, new Set(["g"])).map((row) => row.id)).toEqual(["a", "g", "c"]);

        // c dragged up over b takes its place: first in the group, before b
        expect(projectLayerTreeDrop(rows, "c", "b", 0)).toEqual({ depth: 1, parentId: "g", afterId: null });
        // b dragged left stays where it is but leaves the group, after it
        expect(projectLayerTreeDrop(rows, "b", "b", -1)).toEqual({ depth: 0, parentId: null, afterId: "g" });
        // a dragged down past the group's row lands first in the group: the row below it is the group's
        expect(projectLayerTreeDrop(rows, "a", "g", 0)).toEqual({ depth: 1, parentId: "g", afterId: null });
        // ...unless the group is collapsed, when it goes after the group
        const collapsed = flattenLayerTree(stack.entries, new Set(["g"]));
        expect(projectLayerTreeDrop(collapsed, "a", "g", 1)).toEqual({ depth: 0, parentId: null, afterId: "g" });
    });
});
