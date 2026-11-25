import type { RenderStack, RenderStackEntry } from "@spatialdata/layers";
import type { LayerConfig, RenderStackLayerInputs } from "@spatialdata/vis";
import type { ColorGridSplit } from "./colorGridUtils";

/**
 * A grid of the SpatialData chart by a shapes or labels layer's fill colour: a cell per column
 * the layer colours by, or per category of the first, each drawing its own copy of the layer.
 */
export type LayerFillGrid = { layer: string; split: ColorGridSplit };

type LayerFillGridModeConfig = {
    type?: string;
    grid?: { by: string; split?: string; layer?: string } | null;
};

/** The layer and split of the chart's grid by a layer's fill colour, or null for any other grid. */
export function layerFillGrid(config: LayerFillGridModeConfig): LayerFillGrid | null {
    const { grid } = config;
    if (config.type !== "SpatialDataMdvRegionReact" || grid?.by !== "layer_fill" || !grid.layer) return null;
    return { layer: grid.layer, split: grid.split === "categories" ? "categories" : "columns" };
}

/** The `grid` config of a grid by this layer's fill colour, split this way. */
export function layerFillGridConfig(
    layer: string,
    split: ColorGridSplit,
): { by: "layer_fill"; layer: string; split: ColorGridSplit } {
    return { by: "layer_fill", layer, split };
}

const LAYER_FILL_CHOICE = /^layer_fill:(columns|categories):(.+)$/;

/** A grid settings dropdown's value for a grid by this layer's fill colour. */
export function layerFillChoiceValue(layer: string, split: ColorGridSplit) {
    return `layer_fill:${split}:${layer}`;
}

/** The choices a grid settings dropdown offers for each layer that can be coloured by a column. */
export function layerFillGridChoices(layers: readonly { id: string; name: string }[]) {
    return layers.flatMap(({ id, name }) => [
        { name: `A cell per colour column of ${name}`, value: layerFillChoiceValue(id, "columns") },
        { name: `A cell per colour category of ${name}`, value: layerFillChoiceValue(id, "categories") },
    ]);
}

/** The layer grid a grid settings dropdown's value chooses, or null for another kind of grid. */
export function layerFillGridChoice(value: string): LayerFillGrid | null {
    const match = LAYER_FILL_CHOICE.exec(value);
    if (!match) return null;
    return { layer: match[2], split: match[1] === "categories" ? "categories" : "columns" };
}

/** One cell's copy of the gridded layer: its id (under the cell's view id) and what it draws. */
export type LayerFillCopy = {
    id: string;
    /** The column the copy is coloured by. */
    columnName: string;
    /** Split by categories, the one category of `columnName` whose features the copy draws. */
    category?: string;
};

function isFillColorLayer(
    layer: LayerConfig | undefined,
): layer is Extract<LayerConfig, { type: "shapes" | "labels" }> {
    return layer?.type === "shapes" || layer?.type === "labels";
}

/**
 * The layer inputs with the gridded layer replaced, where it was in the order, by a copy per
 * cell, each coloured by its own column. The copies keep the layer's element, so the renderer
 * loads its geometry once for all of them. Unchanged when the layer is not a shapes or labels
 * layer, or there are no copies.
 */
export function withLayerFillCopies(
    inputs: RenderStackLayerInputs,
    layerId: string,
    copies: readonly LayerFillCopy[],
): RenderStackLayerInputs {
    const layer = inputs.layers[layerId];
    if (!isFillColorLayer(layer) || copies.length === 0) return inputs;
    const { [layerId]: _replaced, ...layers } = inputs.layers;
    const fill = "fillColorByColumn" in layer ? layer.fillColorByColumn : undefined;
    for (const copy of copies) {
        layers[copy.id] = {
            ...layer,
            id: copy.id,
            fillColorByColumn: fill
                ? { ...fill, columnName: copy.columnName }
                : { columnName: copy.columnName, mode: "auto" },
        };
    }
    const layerOrder = inputs.layerOrder.flatMap((id) => (id === layerId ? copies.map((copy) => copy.id) : [id]));
    return { layers, layerOrder };
}

/**
 * The stack in two: what draws below an entry and what draws above it, the entry itself in
 * neither. A group holding the entry is split with it, each half keeping the group's props; the
 * upper half's id is suffixed, as deck needs layer ids to be unique. Null when the entry is not
 * in the stack.
 */
export function splitRenderStackAt(
    stack: RenderStack,
    entryId: string,
): { below: RenderStack; above: RenderStack } | null {
    const split = splitEntries(stack.entries, entryId);
    if (!split) return null;
    return {
        below: { ...stack, entries: split.below },
        above: { ...stack, entries: split.above },
    };
}

function splitEntries(
    entries: readonly RenderStackEntry[],
    entryId: string,
): { below: RenderStackEntry[]; above: RenderStackEntry[] } | null {
    for (const [index, entry] of entries.entries()) {
        const before = entries.slice(0, index);
        const after = entries.slice(index + 1);
        if (entry.id === entryId) return { below: before, above: after };
        if (entry.kind !== "group") continue;
        const inner = splitEntries(entry.children, entryId);
        if (!inner) continue;
        const below = inner.below.length > 0 ? [...before, { ...entry, children: inner.below }] : before;
        const above =
            inner.above.length > 0 ? [{ ...entry, id: `${entry.id}:above`, children: inner.above }, ...after] : after;
        return { below, above };
    }
    return null;
}
