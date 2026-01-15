import { useEffect, useState } from "react";
import { useConfig, useSimplerFilteredIndices } from "../hooks";
import type { TableChartReactConfig } from "../components/TableChartReactWrapper";
import { useDataStore } from "../context";
import { autorun, trace } from "mobx";

// This follows a similar approach for sorting as DataModel.sort()
/**
 * Custom hook to sort and handle externally filtered indices
 * @returns the sorted indices
 */
const useSortedFilteredIndices = () => {
    const config = useConfig<TableChartReactConfig>();
    // const filteredIndices = useReactiveFilteredIndices();
    const filteredIndices = useSimplerFilteredIndices();
    const [sortedFilteredIndices, setSortedFilteredIndices] = useState<Uint32Array>(new Uint32Array(0));
    const dataStore = useDataStore();

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
                    const comparison = decodedData[a].localeCompare(decodedData[b]);
                    return ascending ? comparison : -comparison;
                });
                setSortedFilteredIndices(new Uint32Array(indices));
            } else {
                // All other datatypes
                indices.sort((a, b) => {
                    const valueA = data?.[a];
                    const valueB = data?.[b];

                    // Handle null/undefined (but not 0 or false)
                    const aIsNull = valueA == null;
                    const bIsNull = valueB == null;
                    if (aIsNull && bIsNull) return 0;
                    if (aIsNull) return 1; // null values go to the end
                    if (bIsNull) return -1;

                    // Handle NaN values
                    const aIsNaN = Number.isNaN(valueA);
                    const bIsNaN = Number.isNaN(valueB);
                    if (aIsNaN && bIsNaN) return 0;
                    if (aIsNaN) return 1; // NaN values go to the end
                    if (bIsNaN) return -1;

                    // Normal comparison
                    const comparison = valueA < valueB ? -1 : valueA > valueB ? 1 : 0;
                    return ascending ? comparison : -comparison;
                });
                setSortedFilteredIndices(new Uint32Array(indices));
            }
            return;
        });
        return () => disposer();
    }, [filteredIndices, dataStore, config.sort]);

    return sortedFilteredIndices;
};

export default useSortedFilteredIndices;
