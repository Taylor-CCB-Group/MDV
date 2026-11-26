import type { Layer } from "@deck.gl/core";
import { useMemo } from "react";
import type { DataColumn, DataType } from "@/charts/charts";
import SpatialLayer from "@/webgl/SpatialLayer";
import type { RowColor } from "../colorFields";
import type { UseChartArrayGridResult } from "../hooks/useChartArrayGrid";
import type { ColorGridCells } from "../hooks/useColorGridCells";
import { type CategoryRows, partitionRowsByCategory } from "./colorGridUtils";
import type { DeckCellLabel } from "./DeckCellGrid";
import { subsetPointsLayer } from "./subsetPointsLayer";

const NO_ROWS: CategoryRows = { rows: new Uint32Array(0), positions: new Uint32Array(0) };
const NO_VALUES: readonly string[] = [];
const NO_LAYERS: readonly Layer[] = [];
const NO_HIGHLIGHT: number[] = [];
/** The other categories' points, under a cell's own: darker than filtered-out points. */
const OTHERS_GREY: [number, number, number] = [120, 120, 120];

/** A copy of a points layer under another id, coloured by another column (see `useScatterplotLayer`). */
export type RecolorPoints = (points: Layer, id: string, color: RowColor, column: DataColumn<DataType>) => Layer;

/**
 * A copy of a grey (filtered-out) points layer under another id, leaving out those with no value
 * in another column; null when one grey layer serves every cell (see `useScatterplotLayer`).
 */
export type RegreyPoints = ((grey: Layer, id: string, column: DataColumn<DataType>) => Layer) | null;

/**
 * For a colour grid: the cells in view, each one's points (a copy of `points` with its id prefixed
 * with the cell's view id, so it draws only there) and a label per cell.
 *
 * - Split by columns, a cell's copy draws every point, coloured by its column (`recolorPoints`).
 *   Each copy builds every attribute, so the cost grows with the cells in view. When missing
 *   values are hidden, which differ by column, each cell also draws its own grey points
 *   (`regreyPoints`, under its points), and `ownGrey` says the shared grey layer is to be left out.
 * - Split by categories, a cell's copy draws only its category's rows, and its label says how
 *   many. The rows are split once, by the column's codes, so each copy builds its attributes from
 *   its own rows only. With `greyOthers`, every cell also draws the other categories' points in
 *   grey, under its own: one layer of all the points (`sharedLayers`), drawn in every cell.
 */
export function useColorGridCellLayers(
    grid: UseChartArrayGridResult,
    colorCells: ColorGridCells,
    /** The points layer over all the rows the chart draws; absent while hidden. */
    points: Layer | null | undefined,
    recolorPoints: RecolorPoints,
    /** The grey (filtered-out) points layer, and how to copy it for a cell. */
    grey?: { layer: Layer | null | undefined; regreyPoints: RegreyPoints },
    /** Split by categories: draw the other categories' points in grey in each cell. */
    greyOthers = false,
) {
    const greyLayer = grey?.regreyPoints ? grey.layer : undefined;
    const regreyPoints = greyLayer ? grey?.regreyPoints : undefined;
    const { cells } = colorCells;
    const visibleCells = useMemo(
        () => grid.visibleCellIndices.filter((index) => cells[index] && grid.viewIds[index]),
        [grid.visibleCellIndices, grid.viewIds, cells],
    );
    const codes = colorCells.split === "categories" ? colorCells.codes : undefined;
    const values = colorCells.split === "categories" ? colorCells.values : NO_VALUES;
    const dataVersion = colorCells.split === "categories" ? colorCells.dataVersion : 0;
    const pointRows = points?.props.data;
    // `cells` changes when the column gains categories in place, which `values` would not show, and
    // `dataVersion` when rows move between categories, which none of the arrays would
    // biome-ignore lint/correctness/useExhaustiveDependencies: see above
    const rowsByCategory = useMemo(() => {
        if (!(pointRows instanceof Uint32Array) || !codes) return null;
        return partitionRowsByCategory(pointRows, codes, values);
    }, [pointRows, codes, values, cells, dataVersion]);
    // the highlight is drawn apart from `data`: each cell gets its own category's
    const highlighted = points instanceof SpatialLayer ? points.props.highlightedData : undefined;
    // biome-ignore lint/correctness/useExhaustiveDependencies: see `rowsByCategory`
    const highlightedByCategory = useMemo(() => {
        if (!Array.isArray(highlighted) || highlighted.length === 0 || !codes) return null;
        return partitionRowsByCategory(highlighted, codes, values);
    }, [highlighted, codes, values, cells, dataVersion]);

    const pointLayers = useMemo(() => {
        if (!points) return [];
        return visibleCells.flatMap((index) => {
            const cell = cells[index];
            const viewId = grid.viewIds[index];
            if (!cell || !viewId) return [];
            const id = `${viewId}-${points.id}`;
            if (colorCells.split === "columns") {
                const column = colorCells.columns[index];
                const color = colorCells.colorFunctions[index];
                if (!column || !color) return [];
                const recolored = recolorPoints(points, id, color, column);
                if (!greyLayer || !regreyPoints) return [recolored];
                return [regreyPoints(greyLayer, `${viewId}-${greyLayer.id}`, column), recolored];
            }
            const rows = rowsByCategory?.get(cell.key) ?? NO_ROWS;
            const cellHighlighted = Array.from(highlightedByCategory?.get(cell.key)?.rows ?? []);
            return [subsetPointsLayer(points, rows, cellHighlighted, id)];
        });
    }, [
        points,
        visibleCells,
        cells,
        grid.viewIds,
        colorCells,
        recolorPoints,
        greyLayer,
        regreyPoints,
        rowsByCategory,
        highlightedByCategory,
    ]);
    const ownGrey = colorCells.split === "columns" && !!greyLayer;
    const showOthers = greyOthers && colorCells.split === "categories";
    const sharedLayers = useMemo(() => {
        if (!points || !showOthers) return NO_LAYERS;
        // a points layer's props, which a `Layer` in general does not have
        const props = {
            id: `${points.id}-others`,
            getFillColor: OTHERS_GREY,
            pickable: false,
            ...(points instanceof SpatialLayer && { highlightedData: NO_HIGHLIGHT }),
            updateTriggers: { ...points.props.updateTriggers, getFillColor: OTHERS_GREY },
        };
        return [points.clone(props)];
    }, [points, showOthers]);

    const labels: DeckCellLabel[] = useMemo(
        () =>
            cells.map((cell) => {
                if (!rowsByCategory) return { text: cell.label, color: cell.color };
                const count = rowsByCategory.get(cell.key)?.rows.length ?? 0;
                return { text: `${cell.label} (${count})`, color: cell.color };
            }),
        [cells, rowsByCategory],
    );
    return { visibleCells, pointLayers, labels, ownGrey, sharedLayers };
}
