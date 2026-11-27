import type { SpatialData } from "@spatialdata/core";
import { findRenderStackEntry } from "@spatialdata/layers";
import { action } from "mobx";
import { useEffect, useMemo } from "react";
import { flattenFields } from "@/lib/columnTypeHelpers";
import {
    type ColorGridCell,
    type ColorGridSplit,
    colorGridCells,
    colorGridSplitFor,
    hidesEmptyCells,
} from "../components/colorGridUtils";
import { getDensityGridViewId } from "../components/densityGridUtils";
import {
    fillColumnName,
    type LayerFillCopy,
    layerFillGrid,
    layerFillGridConfig,
} from "../components/layerFillGridUtils";
import type { SpatialDataMdvReactConfig } from "../components/SpatialDataMDVReact";
import { getFieldColor } from "../fieldColorManager";
import { useChartID, useConfig, useDataSources } from "../hooks";
import { mdvFieldSpecsOf } from "../spatialdata/field_spec_projection";
import { resolveAssociatedElementTable, useDataStoreFilterVersion } from "../spatialdata/table_association";

/** The cells of a grid by a layer's fill colour, and the copy of the layer each draws. */
export type LayerFillGridCells = {
    /** The render stack entry the copies stand in for. */
    layerId: string;
    split: ColorGridSplit;
    cells: ColorGridCell[];
    /** Aligned with `cells`. */
    copies: LayerFillCopy[];
    /** Split by categories, each copy's category, by its id (see `useAssociatedShapesLayerInputs`). */
    categoryByLayerId?: Record<string, string>;
};

const NO_VALUES: readonly string[] = [];

/**
 * The cells of the chart's grid by a shapes or labels layer's fill colour, or null when it draws
 * none. Split by columns, a cell per column the layer's colour picker names (its active link's
 * columns as they change); by categories, a cell per category of the first, which needs one
 * category per feature (`text` or `text16`). It keeps the split the first column's type allows. The columns are those of the table annotating the
 * layer's element, which need not be the chart's own datasource.
 */
export function useLayerFillGridCells(spatialData: SpatialData | undefined): LayerFillGridCells | null {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const chartId = useChartID();
    const dataSources = useDataSources();
    const grid = layerFillGrid(config);
    const entry = grid && config.renderStack ? findRenderStackEntry(config.renderStack, grid.layer) : undefined;
    const spatialEntry =
        entry?.kind === "spatial" && (entry.source.elementType === "shapes" || entry.source.elementType === "labels")
            ? entry
            : undefined;
    const elementType = spatialEntry?.source.elementType;
    const elementKey = spatialEntry?.source.elementKey;

    // read here, in an observer's render, so an active link's change of columns re-renders
    const spec = spatialEntry ? mdvFieldSpecsOf(spatialEntry.props)?.fillColorByColumn : undefined;
    const fallbackColumn = spatialEntry ? fillColumnName(spatialEntry.props) : undefined;
    const names = spec ? flattenFields(spec) : fallbackColumn ? [fallbackColumn] : [];
    const namesKey = names.join("\u0000");

    const table = useMemo(
        () =>
            elementType
                ? resolveAssociatedElementTable({ spatialData, elementType, elementKey, dataSources })
                : undefined,
        [spatialData, elementType, elementKey, dataSources],
    );
    const dataStore = table?.status === "resolved" ? table.dataStore : undefined;

    // The split follows the first column's type, when it changes other than by the layer's picker
    // (an active link); none (a link not yet resolved) leaves the grid be.
    const typeSplit = colorGridSplitFor(names[0] && dataStore ? dataStore.columnIndex[names[0]] : undefined);
    const gridLayer = grid?.layer;
    const gridSplit = grid?.split;
    useEffect(() => {
        if (!gridLayer || !typeSplit || typeSplit === gridSplit) return;
        action(() => {
            // a grid chosen since this render (another layer, or none) is left alone
            const current = layerFillGrid(config);
            if (current?.layer !== gridLayer || current.split !== gridSplit) return;
            config.grid = layerFillGridConfig(gridLayer, typeSplit, current);
        })();
    }, [config, gridLayer, gridSplit, typeSplit]);

    // the cells follow the column's type now, before the effect above saves it; the saved split
    // stands while a link has not resolved
    const split = gridLayer ? (typeSplit ?? gridSplit) : undefined;
    const first = split === "categories" ? names[0] : undefined;
    const column = first && dataStore ? dataStore.columnIndex[first] : undefined;
    const categorical = column && (column.datatype === "text" || column.datatype === "text16") ? column : undefined;
    const values: readonly string[] = categorical?.values ?? NO_VALUES;
    // A column gains categories, and changes colours, in place: compare by value.
    const colors =
        categorical && dataStore ? dataStore.getColumnColors(categorical.field, { asArray: true, useValue: true }) : [];
    const valuesKey = values.join("\u0000");
    const colorsKey = JSON.stringify(colors);

    // with `hide_empty`, the categories with no rows left after filtering in the table get no cell
    const hideEmpty = !!categorical && hidesEmptyCells(config.grid);
    const filteredStores = useMemo(() => (hideEmpty && dataStore ? [dataStore] : []), [hideEmpty, dataStore]);
    const filterVersion = useDataStoreFilterVersion(filteredStores);
    // biome-ignore lint/correctness/useExhaustiveDependencies: the filter changes in place (`filterVersion`)
    const presentKey = useMemo(
        () =>
            hideEmpty && dataStore && categorical ? categoriesWithRows(dataStore, categorical.data, values.length) : "",
        [hideEmpty, dataStore, categorical, valuesKey, filterVersion],
    );

    const layerId = gridLayer;
    // biome-ignore lint/correctness/useExhaustiveDependencies: keyed by value, as columns and links change in place
    return useMemo(() => {
        if (!layerId || !split || !spatialEntry) return null;
        const copyId = (cell: ColorGridCell, index: number) =>
            `${getDensityGridViewId(chartId, cell.key, index)}-${layerId}`;
        if (split === "categories") {
            if (!first) return null;
            const all = colorGridCells([...values], colors);
            const cells = presentKey ? all.filter((_cell, index) => presentKey[index] === "1") : all;
            const copies = cells.map((cell, index) => ({
                id: copyId(cell, index),
                columnName: first,
                category: cell.key,
            }));
            const categoryByLayerId = Object.fromEntries(copies.map((copy) => [copy.id, copy.category]));
            return { layerId, split, cells, copies, categoryByLayerId };
        }
        const unique = Array.from(new Set(names));
        const cells = unique.map((name) => ({
            key: name,
            label: dataStore?.columnIndex[name]?.name || name,
            color: getFieldColor(name),
        }));
        const copies = cells.map((cell, index) => ({ id: copyId(cell, index), columnName: cell.key }));
        return { layerId, split, cells, copies };
    }, [layerId, split, chartId, spatialEntry?.id, first, valuesKey, colorsKey, namesKey, dataStore, presentKey]);
}

/**
 * Which categories have a row that passes the table's filters, as a "1" or "0" per category (a
 * string, so an unchanged answer compares equal). Empty while the column's data has not loaded.
 */
function categoriesWithRows(
    dataStore: { size: number; filterArray: ArrayLike<number> },
    codes: ArrayLike<number> | undefined,
    categoryCount: number,
) {
    if (!codes) return "";
    const present = new Uint8Array(categoryCount);
    for (let row = 0; row < dataStore.size; row++) {
        if (dataStore.filterArray[row] !== 0) continue;
        const code = codes[row];
        if (code < categoryCount) present[code] = 1;
    }
    return present.join("");
}
