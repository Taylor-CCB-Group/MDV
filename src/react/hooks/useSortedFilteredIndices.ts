import { useEffect, useState } from "react";
import { useConfig, useSimplerFilteredIndices } from "../hooks";
import type { TableChartReactConfig } from "../components/TableChartReactWrapper";
import { useDataStore } from "../context";
import { autorun } from "mobx";

/**
 * Hook that sorts the filtered indices based on the config.sort
 * 
 * - Follows the logic of sorting in DataModel.sort()
 * - Uses Mobx autorun to react to config.sort changes
 * - Updates the indices when any of these change: filteredIndices,
 * dataStore or config.sort
 * 
 * For unique columns, we decode the data and sort it
 * For all other columns we directly sort it, we put the null and 
 * NaN values at the end
 * 
 * Returns a new Uint32Array of sorted indices
 */
const useSortedFilteredIndices = () => {
    const config = useConfig<TableChartReactConfig>();
    // const filteredIndices = useReactiveFilteredIndices();
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
                        `Column ${columnId} of type 'unique' has invalid or missing stringLength: ${length}.`,
                    );
                    setSortedFilteredIndices(indices);
                    return;
                }

                // Decode the data
                for (let i = 0; i < indices.length; i++) {
                    const dataIndex = indices[i];

                    if (!data || dataIndex * length + length > data.length) {
                        console.error(`Index out of bounds for the column ${columnId}, skipping.`);
                        decodedData[dataIndex] = "";
                        continue;
                    }

                    decodedData[dataIndex] = textDecoder.decode(
                        data?.slice?.(dataIndex * length, dataIndex * length + length),
                    );
                }

                // Sort the indices by comparing decoded strings
                indices.sort((a, b) => {
                    const strA = decodedData[a];
                    const strB = decodedData[b];

                    // Handle null/undefined and empty strings
                    const aIsEmpty = !strA || strA.trim() === "";
                    const bIsEmpty = !strB || strB.trim() === "";

                    if (aIsEmpty && bIsEmpty) return 0;
                    if (aIsEmpty) return 1;
                    if (bIsEmpty) return -1;

                    const comparison = strA.localeCompare(strB);
                    return ascending ? comparison : -comparison;
                });
                setSortedFilteredIndices(new Uint32Array(indices));
            } else {
                // All other datatypes
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
            return;
        });
        return () => disposer();
    }, [filteredIndices, dataStore]);

    return sortedFilteredIndices;
};

export default useSortedFilteredIndices;
