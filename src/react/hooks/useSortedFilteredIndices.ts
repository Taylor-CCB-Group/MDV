import { useEffect, useState } from "react";
import { useConfig, useSimplerFilteredIndices } from "../hooks";
import type { TableChartReactConfig } from "../components/TableChartReactWrapper";
import { useDataStore } from "../context";
import { autorun } from "mobx";

/**
 * 
 * @param indices - array of indices to be sorted
 * @param decodedData - map of decoded data
 * @param ascending - order of sort
 * @param isEmpty - function to check if the value is empty or not
 * 
 * Uses Collator to naturally sort the string, this is useful when 
 * there is a numeric part in the string
 */
function sortDecodedStrings(
    indices: Uint32Array,
    decodedData: Record<number, string>,
    ascending: boolean,
    isEmpty: (str: string) => boolean,
): void {
    // Cursor's recommendation to naturally sort the string when there is a numeric part in it
    const collator = new Intl.Collator(undefined, { numeric: true });
    indices.sort((a, b) => {
        const strA = decodedData[a];
        const strB = decodedData[b];
        const aIsEmpty = isEmpty(strA);
        const bIsEmpty = isEmpty(strB);
        if (aIsEmpty && bIsEmpty) return 0;
        if (aIsEmpty) return 1;
        if (bIsEmpty) return -1;
        const comparison = collator.compare(strA, strB);
        return ascending ? comparison : -comparison;
    });
}
/**
 * Hook that sorts the filtered indices based on the config.sort
 *
 * - Follows the logic of sorting in DataModel.sort()
 * - Uses Mobx autorun to react to config.sort changes
 * - Updates the indices when any of these change: filteredIndices,
 * dataStore or config.sort
 *
 * For unique, multitext, text and text16 columns, we decode the data and sort by string.
 * For all other columns we directly sort by raw value.
 * The null, NaN and empty values are put at the end.
 * 
 * ! Multitext column values are sorted based on the displayed values currently
 * Ex: ["B, B", "C, C", "A, A"] -> Sorted: ["A, A", "B, B", "C, C"]
 *
 * Returns a new Uint32Array of sorted indices
 */
const useSortedFilteredIndices = () => {
    const config = useConfig<TableChartReactConfig>();
    const filteredIndices = useSimplerFilteredIndices();
    const [sortedFilteredIndices, setSortedFilteredIndices] = useState<Uint32Array>(new Uint32Array(0));
    const dataStore = useDataStore();

    // If we include the config.sort in the useEffect's dependencies, then it will cause unnecessary rerenders when config.sort changes everytime
    // biome-ignore lint/correctness/useExhaustiveDependencies: the config.sort is tracked by mobx autorun
    useEffect(() => {
        const disposer = autorun(() => {
            // Create a copy of filtered indices
            const indices = new Uint32Array(filteredIndices);

            const sortConfig = config.sort;
            if (!sortConfig) {
                setSortedFilteredIndices(indices);
                return;
            }

            const { columnId, ascending } = sortConfig;

            // Handle index column
            if (columnId === "__index__") {
                indices.sort();
                if (!ascending) {
                    indices.reverse();
                }
                setSortedFilteredIndices(indices);
                return;
            }

            // Get column info
            const colInfo = dataStore.columnIndex[columnId];
            if (!colInfo) {
                setSortedFilteredIndices(indices);
                return;
            }

            // Get raw column data
            const data = dataStore.getRawColumn(columnId);

            if (colInfo.datatype === "unique") {
                // Text sorting for unique columns
                const decodedData: Record<number, string> = {};
                const textDecoder = new TextDecoder();
                const length = colInfo.stringLength;

                if (!length || typeof length !== "number" || length <= 0) {
                    console.error(
                        `Column ${columnId} of type ${colInfo.datatype} has invalid or missing stringLength: ${length}.`,
                    );
                    setSortedFilteredIndices(indices);
                    return;
                }

                // Decode the data
                for (let i = 0; i < indices.length; i++) {
                    // Get the data index from the row indices array
                    const dataIndex = indices[i];

                    if (!data || dataIndex * length + length > data.length) {
                        console.error(`Index out of bounds for the column ${columnId}, skipping.`);
                        decodedData[dataIndex] = "";
                        continue;
                    }

                    // decode the data and assign it to the map with data index as key
                    decodedData[dataIndex] = textDecoder.decode(
                        data?.slice?.(dataIndex * length, dataIndex * length + length),
                    );
                }

                // Sort the indices by comparing decoded strings
                sortDecodedStrings(indices, decodedData, ascending, (s) => !s || s.trim() === "");
                setSortedFilteredIndices(new Uint32Array(indices));

            } else if (colInfo.datatype === "multitext") {
                const decodedData: Record<number, string> = {};
                const length = colInfo.stringLength;

                if (!length || typeof length !== "number" || length <= 0) {
                    console.error(
                        `Column ${columnId} of type ${colInfo.datatype} has invalid or missing stringLength: ${length}.`,
                    );
                    setSortedFilteredIndices(indices);
                    return;
                }

                for (let i = 0; i < indices.length; i++) {
                    // Get the data index from the row indices array
                    const dataIndex = indices[i];
                    if (!data || dataIndex * length + length > data.length) {
                        console.error(`Index out of bounds for the column ${columnId}, skipping.`);
                        decodedData[dataIndex] = "";
                        continue;
                    }

                    // Get the display value of the row by joining the values by the separator - ', '
                    const rowValue = Array.from(data?.slice?.(dataIndex * length, dataIndex * length + length))
                        .filter((x) => x !== 65535)
                        .map((x) => colInfo.values[x] ?? "")
                        .join(", ");
                    decodedData[dataIndex] = rowValue;
                }

                // Sort decoded displayed strings
                sortDecodedStrings(indices, decodedData, ascending, (s) => !s || s.trim() === "" || s.trim() === "N/A");
                setSortedFilteredIndices(new Uint32Array(indices));

            } else if (colInfo.datatype === "text" || colInfo.datatype === "text16") {

                const values = colInfo.values;

                if (!values || !Array.isArray(values)) {
                    setSortedFilteredIndices(indices);
                    return;
                }

                const decodedData: Record<number, string> = {};

                for (let i = 0; i < indices.length; i++) {
                    // Get the data index (data array index)
                    const dataIndex = indices[i];
                    // Get the value index (values array index)
                    const valueIndex = data?.[dataIndex];
                    
                    const isValidValueIndex = valueIndex != null && valueIndex >= 0 && valueIndex < values.length
                    decodedData[dataIndex] =
                        isValidValueIndex
                            ? values[valueIndex] ?? ""
                            : "";
                }

                sortDecodedStrings(indices, decodedData, ascending, (s) => !s || s.trim() === "");
                setSortedFilteredIndices(new Uint32Array(indices));

            } else {
                // All other datatypes (numeric, etc.)
                indices.sort((a, b) => {
                    const valueA = data?.[a];
                    const valueB = data?.[b];

                    // Handle null/undefined and NaN values (but not 0 or false)
                    const aIsEmpty = valueA == null || Number.isNaN(valueA);
                    const bIsEmpty = valueB == null || Number.isNaN(valueB);
                    if (aIsEmpty && bIsEmpty) return 0;
                    if (aIsEmpty) return 1;
                    if (bIsEmpty) return -1;

                    // Normal comparison
                    const comparison = valueA < valueB ? -1 : valueA > valueB ? 1 : 0;
                    return ascending ? comparison : -comparison;
                });
                setSortedFilteredIndices(new Uint32Array(indices));
            }
        });
        return () => disposer();
    }, [filteredIndices, dataStore]);

    return sortedFilteredIndices;
};

export default useSortedFilteredIndices;
