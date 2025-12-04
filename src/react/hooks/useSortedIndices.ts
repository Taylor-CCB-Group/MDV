import { useEffect, useState } from "react";
import { useConfig, useSimplerFilteredIndices } from "../hooks";
import type { TableChartReactConfig } from "../components/TableChartReactWrapper";
import { useDataStore } from "../context";
import { autorun, trace } from "mobx";

// This follows a similar approach for sorting as DataModel.sort()
const useSortedIndices = () => {
    const config = useConfig<TableChartReactConfig>();
    const filteredIndices = useSimplerFilteredIndices();
    const [sortedIndices, setSortedIndices] = useState<Uint32Array>(new Uint32Array(0));
    const dataStore = useDataStore();

    useEffect(() => {
        const disposer = autorun(() => {
            // Create a copy of filtered indices
            const indices = new Uint32Array(filteredIndices);

            const sortConfig = config.sort;
            if (!sortConfig) {
                setSortedIndices(indices);
                return;
            }

            const { columnId, ascending } = sortConfig;

            // Handle index column
            if (columnId === "__index__") {
                indices.sort();
                if (!ascending) {
                    indices.reverse();
                }
                setSortedIndices(indices);
                return;
            }

            // Get column info
            const colInfo = dataStore.columnIndex[columnId];
            if (!colInfo) {
                setSortedIndices(indices);
                return;
            }

            // Get raw column data
            const data = dataStore.getRawColumn(columnId);

            if (colInfo.datatype === "unique") {
                // Text sorting for unique columns
                const decodedData: Record<number, string> = {};
                const textDecoder = new TextDecoder();
                const length = colInfo.stringLength;

                // Decode the data
                for (let i = 0; i < indices.length; i++) {
                    const dataIndex = indices[i];
                    decodedData[dataIndex] = textDecoder.decode(
                        data?.slice?.(dataIndex * length, dataIndex * length + length)
                    );
                }

                // Sort the indices by comparing decoded strings
                const sorted = Array.from(indices).sort((a, b) => {
                    const comparison = decodedData[a].localeCompare(decodedData[b]);
                    return ascending ? comparison : -comparison;
                });
                setSortedIndices(new Uint32Array(sorted));
            } else {
                const sorted = Array.from(indices).sort((a, b) => {
                    const valueA = data?.[a];
                    const valueB = data?.[b];

                    if (!valueA || !valueB) {
                        return 1;
                    }

                    // Handle NaN values
                    const aIsNaN = Number.isNaN(valueA);
                    const bIsNaN = Number.isNaN(valueB);
                    if (aIsNaN && bIsNaN) return 0;
                    if (aIsNaN) return 1;
                    if (bIsNaN) return -1;

                    const comparison = valueA < valueB ? -1 : valueA > valueB ? 1 : 0;
                    return ascending ? comparison : -comparison;
                });
                setSortedIndices(new Uint32Array(sorted));
            }
            return;
        });
        return () => disposer();
    }, [filteredIndices, dataStore, config.sort]);

    return sortedIndices;

};

export default useSortedIndices;