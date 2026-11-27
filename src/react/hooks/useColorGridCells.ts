import { useEffect, useId, useMemo, useState } from "react";
import type { DataColumn, DataType } from "@/charts/charts";
import type DataStore from "@/datastore/DataStore";
import { type ColorFieldsConfig, type RowColor, useColorFields } from "../colorFields";
import { type ColorGridCell, colorGridCells, colorGridSplit } from "../components/colorGridUtils";
import { useDataStore } from "../context";
import { getFieldColor } from "../fieldColorManager";
import { useConfig } from "../hooks";

/** The cells of a colour grid, by how it splits the colour. */
export type ColorGridCells =
    | {
          /** A cell per category of the first colour column, each drawing that category's points. */
          split: "categories";
          cells: ColorGridCell[];
          /** The column's data, a code per row indexing `values`: which cell a row is drawn in. */
          codes: ArrayLike<number> | undefined;
          values: readonly string[];
          /** Advances when data is changed in place: rows can move category with the same arrays. */
          dataVersion: number;
      }
    | {
          /** A cell per colour column, each drawing every point coloured by its column. */
          split: "columns";
          cells: ColorGridCell[];
          /** Aligned with `cells`. */
          columns: DataColumn<DataType>[];
          /** Aligned with `cells`. */
          colorFunctions: RowColor[];
      };

const NO_VALUES: readonly string[] = [];

/** A version that advances when the datastore's data is changed in place (`data_changed`). */
function useDataChangedVersion(dataStore: DataStore, active: boolean) {
    const [version, setVersion] = useState(0);
    const instanceId = useId();
    useEffect(() => {
        if (!active) return;
        const listenerId = `color-grid-cells-${instanceId}`;
        dataStore.addListener(listenerId, (type: string) => {
            if (type === "data_changed") setVersion((current) => current + 1);
        });
        return () => dataStore.removeListener(listenerId);
    }, [dataStore, active, instanceId]);
    return version;
}

/**
 * The cells of the chart's colour grid, or null when it draws none. Split by columns, a cell per
 * loaded colour column; by categories, a cell per category of the first, which needs one category
 * per row (`text` or `text16`). A cell is made only once its column has loaded.
 */
export function useColorGridCells(): ColorGridCells | null {
    const config = useConfig<ColorFieldsConfig>();
    const dataStore = useDataStore();
    const split = colorGridSplit(config);
    const { columns, colorFunctions } = useColorFields();
    const dataVersion = useDataChangedVersion(dataStore, split === "categories");

    const first = split === "categories" ? columns[0] : undefined;
    const categorical = first?.data && (first.datatype === "text" || first.datatype === "text16") ? first : null;
    const values = categorical?.values ?? NO_VALUES;
    // A column gains categories, and changes colours, in place: compare by value.
    const colors = categorical ? dataStore.getColumnColors(categorical.field, { asArray: true, useValue: true }) : [];
    const valuesKey = values.join("\u0000");
    const colorsKey = JSON.stringify(colors);
    // biome-ignore lint/correctness/useExhaustiveDependencies: keyed by value, as the column changes in place
    const categoryCells = useMemo(() => colorGridCells([...values], colors), [valuesKey, colorsKey]);

    const columnCells = useMemo(
        () =>
            split === "columns"
                ? columns.map((column) => ({
                      key: column.field,
                      label: column.name || column.field,
                      color: getFieldColor(column.field),
                  }))
                : [],
        [split, columns],
    );

    if (split === "categories") {
        return { split, cells: categoryCells, codes: categorical?.data, values, dataVersion };
    }
    if (split === "columns") return { split, cells: columnCells, columns, colorFunctions };
    return null;
}
