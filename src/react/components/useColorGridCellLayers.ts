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

/** A copy of a points layer under another id, coloured by another column (see `useScatterplotLayer`). */
export type RecolorPoints = (points: Layer, id: string, color: RowColor, column: DataColumn<DataType>) => Layer;

/**
 * For a colour grid: the cells in view, each one's points (a copy of `points` with its id prefixed
 * with the cell's view id, so it draws only there) and a label per cell.
 *
 * - Split by columns, a cell's copy draws every point, coloured by its column (`recolorPoints`).
 *   Each copy builds every attribute, so the cost grows with the cells in view.
 * - Split by categories, a cell's copy draws only its category's rows, and its label says how
 *   many. The rows are split once, by the column's codes, so each copy builds its attributes from
 *   its own rows only.
 */
export function useColorGridCellLayers(
    grid: UseChartArrayGridResult,
    colorCells: ColorGridCells,
    /** The points layer over all the rows the chart draws; absent while hidden. */
    points: Layer | null | undefined,
    recolorPoints: RecolorPoints,
) {
    const { cells } = colorCells;
    const visibleCells = useMemo(
        () => grid.visibleCellIndices.filter((index) => cells[index] && grid.viewIds[index]),
        [grid.visibleCellIndices, grid.viewIds, cells],
    );
    const codes = colorCells.split === "categories" ? colorCells.codes : undefined;
    const values = colorCells.split === "categories" ? colorCells.values : NO_VALUES;
    const pointRows = points?.props.data;
    // `cells` changes when the column gains categories in place, which `values` would not show
    // biome-ignore lint/correctness/useExhaustiveDependencies: see above
    const rowsByCategory = useMemo(() => {
        if (!(pointRows instanceof Uint32Array) || !codes) return null;
        return partitionRowsByCategory(pointRows, codes, values);
    }, [pointRows, codes, values, cells]);
    // the highlight is drawn apart from `data`: each cell gets its own category's
    const highlighted = points instanceof SpatialLayer ? points.props.highlightedData : undefined;
    // biome-ignore lint/correctness/useExhaustiveDependencies: see `rowsByCategory`
    const highlightedByCategory = useMemo(() => {
        if (!Array.isArray(highlighted) || highlighted.length === 0 || !codes) return null;
        return partitionRowsByCategory(highlighted, codes, values);
    }, [highlighted, codes, values, cells]);

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
                return column && color ? [recolorPoints(points, id, color, column)] : [];
            }
            const rows = rowsByCategory?.get(cell.key) ?? NO_ROWS;
            const cellHighlighted = Array.from(highlightedByCategory?.get(cell.key)?.rows ?? []);
            return [subsetPointsLayer(points, rows, cellHighlighted, id)];
        });
    }, [points, visibleCells, cells, grid.viewIds, colorCells, recolorPoints, rowsByCategory, highlightedByCategory]);

    const labels: DeckCellLabel[] = useMemo(
        () =>
            cells.map((cell) => {
                if (!rowsByCategory) return { text: cell.label, color: cell.color };
                const count = rowsByCategory.get(cell.key)?.rows.length ?? 0;
                return { text: `${cell.label} (${count})`, color: cell.color };
            }),
        [cells, rowsByCategory],
    );
    return { visibleCells, pointLayers, labels };
}
