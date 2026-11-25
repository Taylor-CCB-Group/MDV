import type { RenderStack, RenderStackEntry } from "@spatialdata/layers";
import type { LayerConfig, RenderStackLayerInputs } from "@spatialdata/vis";
import { describe, expect, test } from "vitest";
import { SpatialDataGridConfigSchema } from "@/charts/schemas/ChartConfigSchema";
import { rowsInCategory } from "@/react/spatialdata/table_association";
import { gridChoiceValue } from "./colorGridUtils";
import {
    layerFillGrid,
    layerFillGridChoice,
    layerFillGridChoices,
    layerFillGridConfig,
    splitRenderStackAt,
    withLayerFillCopies,
} from "./layerFillGridUtils";

const SPATIAL = "SpatialDataMdvRegionReact";

describe("layerFillGrid", () => {
    test("reads the layer and split of a SpatialData chart's grid by a layer's fill", () => {
        expect(layerFillGrid({ type: SPATIAL, grid: { by: "layer_fill", layer: "cells" } })).toEqual({
            layer: "cells",
            split: "columns",
        });
        expect(
            layerFillGrid({ type: SPATIAL, grid: { by: "layer_fill", layer: "cells", split: "categories" } }),
        ).toEqual({ layer: "cells", split: "categories" });
    });

    test("is null for another grid, another chart, or no layer", () => {
        expect(layerFillGrid({ type: SPATIAL, grid: { by: "color_fields" } })).toBeNull();
        expect(layerFillGrid({ type: "DeckScatter", grid: { by: "layer_fill", layer: "cells" } })).toBeNull();
        expect(layerFillGrid({ type: SPATIAL, grid: { by: "layer_fill" } })).toBeNull();
        expect(layerFillGrid({ type: SPATIAL, grid: null })).toBeNull();
    });
});

describe("layer grid choices", () => {
    test("round trip through the grid settings dropdown", () => {
        const id = "spatialdata-shapes-a-#x1";
        const choices = layerFillGridChoices([
            { id, name: "shapes: a", split: "columns" },
            { id, name: "shapes: a", split: "categories" },
        ]);
        expect(choices.map((choice) => choice.name)).toEqual([
            "A cell per colour column of shapes: a",
            "A cell per colour category of shapes: a",
        ]);
        for (const choice of choices) {
            const grid = layerFillGridChoice(choice.value);
            expect(grid).not.toBeNull();
            if (!grid) continue;
            expect(gridChoiceValue({ grid: layerFillGridConfig(grid.layer, grid.split) })).toBe(choice.value);
        }
        expect(layerFillGridChoice("color_fields:columns")).toBeNull();
    });

    test("the schema keeps the layer", () => {
        const grid = { by: "layer_fill", layer: "cells", split: "categories" };
        expect(SpatialDataGridConfigSchema.safeParse(grid).data).toEqual(grid);
    });
});

describe("withLayerFillCopies", () => {
    const shapes: LayerConfig = {
        type: "shapes",
        id: "cells",
        elementKey: "cell_circles",
        visible: true,
        opacity: 1,
        fillColorByColumn: { columnName: "a", mode: "auto" },
    };
    const image: LayerConfig = { type: "image", id: "he", elementKey: "he", visible: true, opacity: 1 };
    const inputs: RenderStackLayerInputs = { layers: { he: image, cells: shapes }, layerOrder: ["he", "cells"] };

    test("puts a copy per cell where the layer was, each coloured by its own column", () => {
        const next = withLayerFillCopies(inputs, "cells", [
            { id: "v0-cells", columnName: "a" },
            { id: "v1-cells", columnName: "b" },
        ]);
        expect(next.layerOrder).toEqual(["he", "v0-cells", "v1-cells"]);
        expect(next.layers.cells).toBeUndefined();
        expect(next.layers["v1-cells"]).toMatchObject({
            id: "v1-cells",
            elementKey: "cell_circles",
            fillColorByColumn: { columnName: "b", mode: "auto" },
        });
        expect(next.layers.he).toBe(image);
    });

    test("leaves the inputs alone for a layer that cannot be coloured by a column, or no copies", () => {
        expect(withLayerFillCopies(inputs, "he", [{ id: "v0-he", columnName: "a" }])).toBe(inputs);
        expect(withLayerFillCopies(inputs, "cells", [])).toBe(inputs);
    });
});

describe("splitRenderStackAt", () => {
    const leaf = (id: string): RenderStackEntry => ({
        kind: "host",
        id,
        visible: true,
        props: {},
        source: { hostLayerId: id },
    });
    const stack: RenderStack = {
        schemaVersion: 2,
        entries: [
            leaf("image"),
            {
                kind: "group",
                id: "group",
                visible: true,
                props: { opacity: 0.5 },
                children: [leaf("under"), leaf("cells"), leaf("over")],
            },
            leaf("points"),
        ],
    };

    test("splits the stack, and any group holding the entry, either side of it", () => {
        const split = splitRenderStackAt(stack, "cells");
        expect(split?.below.entries.map((entry) => entry.id)).toEqual(["image", "group"]);
        expect(split?.above.entries.map((entry) => entry.id)).toEqual(["group:above", "points"]);
        const [, lowerGroup] = split?.below.entries ?? [];
        const [upperGroup] = split?.above.entries ?? [];
        expect(lowerGroup?.kind === "group" && lowerGroup.children.map((entry) => entry.id)).toEqual(["under"]);
        expect(upperGroup?.kind === "group" && upperGroup.children.map((entry) => entry.id)).toEqual(["over"]);
        expect(upperGroup?.props).toEqual({ opacity: 0.5 });
    });

    test("is null for an entry not in the stack", () => {
        expect(splitRenderStackAt(stack, "missing")).toBeNull();
    });
});

describe("rowsInCategory", () => {
    test("keeps the rows of one category, in order", () => {
        const codes = [0, 1, 0, 2, 1];
        expect([...rowsInCategory(new Uint32Array([4, 0, 2, 1]), codes, 1)]).toEqual([4, 1]);
        expect([...rowsInCategory([0, 2], codes, 2)]).toEqual([]);
    });
});
