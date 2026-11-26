import type { RenderStack, RenderStackEntry } from "@spatialdata/layers";
import type { LayerConfig, RenderStackLayerInputs } from "@spatialdata/vis";
import { deckHostLayerId } from "../spatialdata/host_overlay_ids";
import { type ColorGridSplit, colorGridSplit } from "./colorGridUtils";
import { isDensityGridShown } from "./densityGridUtils";

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

/**
 * The layer list entry the chart's grid is made from, for marking it there: a layer whose fill
 * colour it splits, the MDV points whose colour it splits (`color_fields`), or the density
 * overlay (a cell per density field). Null when there is no such grid, as for a region grid.
 */
export type GriddedEntry = {
    entryId: string;
    /** Whether the cells split the entry's colour, rather than the entry itself. */
    byColour: boolean;
    /** What each cell is, as words: "a cell per category". */
    cells: string;
};

const CELLS_BY_SPLIT: Record<ColorGridSplit, string> = {
    columns: "a cell per colour column",
    categories: "a cell per category",
};

export function griddedEntry(config: LayerFillGridModeConfig & { dimension?: "2d" | "3d" }): GriddedEntry | null {
    if (config.type !== "SpatialDataMdvRegionReact") return null;
    const layerGrid = layerFillGrid(config);
    if (layerGrid) return { entryId: layerGrid.layer, byColour: true, cells: CELLS_BY_SPLIT[layerGrid.split] };
    const split = colorGridSplit(config);
    if (split) return { entryId: deckHostLayerId("scatter"), byColour: true, cells: CELLS_BY_SPLIT[split] };
    if (isDensityGridShown(config)) {
        return { entryId: deckHostLayerId("density"), byColour: false, cells: "a cell per density field" };
    }
    return null;
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

/** The choices a grid settings dropdown offers for layers, a choice per layer and split. */
export function layerFillGridChoices(layers: readonly { id: string; name: string; split: ColorGridSplit }[]) {
    return layers.map(({ id, name, split }) => ({
        name: split === "columns" ? `A cell per colour column of ${name}` : `A cell per colour category of ${name}`,
        value: layerFillChoiceValue(id, split),
    }));
}

/** The column a shapes or labels layer's props colour it by (the first, for several), if any. */
export function fillColumnName(props: object): string | undefined {
    const fill = "fillColorByColumn" in props ? props.fillColorByColumn : undefined;
    if (!fill || typeof fill !== "object" || !("columnName" in fill)) return undefined;
    const { columnName } = fill;
    return typeof columnName === "string" && columnName ? columnName : undefined;
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
        // sliced only once found, not for every entry passed over
        if (entry.id === entryId) return { below: entries.slice(0, index), above: entries.slice(index + 1) };
        if (entry.kind !== "group") continue;
        const inner = splitEntries(entry.children, entryId);
        if (!inner) continue;
        const before = entries.slice(0, index);
        const after = entries.slice(index + 1);
        const below = inner.below.length > 0 ? [...before, { ...entry, children: inner.below }] : before;
        const above =
            inner.above.length > 0 ? [{ ...entry, id: `${entry.id}:above`, children: inner.above }, ...after] : after;
        return { below, above };
    }
    return null;
}
