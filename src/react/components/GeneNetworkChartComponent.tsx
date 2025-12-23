import { useMemo, useRef, useEffect, useLayoutEffect, useState, type MouseEvent } from "react";
import { observer } from "mobx-react-lite";
import { useChart } from "../context";
import { useParamColumns } from "../hooks";
import { useHighlightedIndices, useHighlightRows } from "../selectionHooks";
import { useSimplerFilteredIndices } from "../hooks";
import { isColumnLoaded } from "@/lib/columnTypeHelpers";
import { GeneNetworkInfoComponent } from "./GeneNetworkInfoComponent";
import type { GeneNetworkConfig } from "./GeneNetworkChart";
import { Chip } from "@mui/material";

/**
 * Get text value from a column at a given index.
 * Uses column.getValue which handles all column types (text, text16, unique, etc.)
 */
function getTextValue(column: any, index: number): string {
    if (!isColumnLoaded(column)) {
        throw new Error("Column not loaded");
    }
    
    // getValue handles all column types correctly
    const value = column.getValue(index);
    return String(value);
}

export const GeneNetworkChartComponent = observer(() => {
    const chart = useChart<GeneNetworkConfig>();
    const columns = useParamColumns();
    const highlightedIndices = useHighlightedIndices();
    const filteredIndices = useSimplerFilteredIndices();
    const highlightRows = useHighlightRows();
    const scrollContainerRef = useRef<HTMLDivElement>(null);
    const geneRefs = useRef<Map<string, HTMLDivElement>>(new Map());
    const [lastClickedIndex, setLastClickedIndex] = useState<number | null>(null);

    if (columns.length === 0) {
        return <div>No column selected. Please configure the chart.</div>;
    }

    const column = columns[0];
    if (!isColumnLoaded(column)) {
        return <div>Loading column data...</div>;
    }

    const mode = chart.config.mode || "filtered";
    const maxGenes = chart.config.maxGenes ?? 100;

    // Get gene IDs based on display mode
    const { geneIds, geneHighlightCounts } = useMemo(() => {
        const highlightCounts = new Map<string, number>();
        
        // First, count how many highlighted rows each gene has
        const highlightedGeneIds = new Set<string>();
        for (let i = 0; i < highlightedIndices.length; i++) {
            const index = highlightedIndices[i];
            try {
                const geneId = getTextValue(column, index);
                if (geneId) {
                    highlightedGeneIds.add(geneId);
                    highlightCounts.set(geneId, (highlightCounts.get(geneId) || 0) + 1);
                }
            } catch (e) {
                // Skip invalid indices
            }
        }

        if (mode === "observableFields") {
            // ObservableFields mode: show highlighted if any, otherwise filtered
            if (highlightedIndices.length > 0) {
                return {
                    geneIds: Array.from(highlightedGeneIds),
                    geneHighlightCounts: highlightCounts,
                };
            }
            // Fall through to filtered logic
        }

        // Filtered mode: get unique gene IDs from filtered indices
        const geneIdSet = new Set<string>();
        try {
            for (let i = 0; i < filteredIndices.length; i++) {
                const index = filteredIndices[i];
                const geneId = getTextValue(column, index);
                if (geneId) {
                    geneIdSet.add(geneId);
                }
            }
        } catch (e) {
            console.error("Error getting filtered gene IDs:", e);
        }
        
        return {
            geneIds: Array.from(geneIdSet),
            geneHighlightCounts: highlightCounts,
        };
    }, [mode, highlightedIndices, filteredIndices, column]);

    // Apply cap
    const visibleGeneIds = geneIds.slice(0, maxGenes);
    const hasMore = geneIds.length > maxGenes;

    // Auto-scroll to highlighted genes
    useLayoutEffect(() => {
        if (highlightedIndices.length === 0 || visibleGeneIds.length === 0) return;
        if (!scrollContainerRef.current) return;

        // Find which genes correspond to highlighted rows
        const highlightedGeneIds = new Set<string>();
        for (let i = 0; i < highlightedIndices.length; i++) {
            const index = highlightedIndices[i];
            try {
                const geneId = getTextValue(column, index);
                if (geneId && visibleGeneIds.includes(geneId)) {
                    highlightedGeneIds.add(geneId);
                }
            } catch (e) {
                // Skip invalid indices
            }
        }

        // Scroll to the first highlighted gene that's not visible
        for (const geneId of highlightedGeneIds) {
            const geneElement = geneRefs.current.get(geneId);
            if (geneElement) {
                const container = scrollContainerRef.current;
                const containerRect = container.getBoundingClientRect();
                const elementRect = geneElement.getBoundingClientRect();

                // Check if element is not fully visible
                const isVisible =
                    elementRect.top >= containerRect.top &&
                    elementRect.bottom <= containerRect.bottom;

                if (!isVisible) {
                    geneElement.scrollIntoView({
                        behavior: "smooth",
                        block: "nearest",
                    });
                    break; // Only scroll to first one
                }
            }
        }
    }, [highlightedIndices, visibleGeneIds, column]);

    const setGeneRef = (geneId: string, element: HTMLDivElement | null) => {
        if (element) {
            geneRefs.current.set(geneId, element);
        } else {
            geneRefs.current.delete(geneId);
        }
    };

    const getIndicesForGeneId = (geneId: string): number[] => {
        const indices: number[] = [];
        for (let i = 0; i < filteredIndices.length; i++) {
            const index = filteredIndices[i];
            try {
                const value = getTextValue(column, index);
                if (value === geneId) {
                    indices.push(index);
                }
            } catch {
                // Skip invalid indices
            }
        }
        return indices;
    };

    const handleGeneClick = (geneId: string, event: MouseEvent<HTMLDivElement>) => {
        event.preventDefault();
        event.stopPropagation();

        const geneRowIndices = getIndicesForGeneId(geneId);
        if (geneRowIndices.length === 0) return;

        const isToggle = event.ctrlKey || event.metaKey;
        const isRange = event.shiftKey;

        // Shift-click: select a continuous range of row indices
        if (isRange) {
            const currentIndex = Math.min(...geneRowIndices);

            if (lastClickedIndex == null) {
                highlightRows(geneRowIndices);
                setLastClickedIndex(currentIndex);
                return;
            }

            const start = Math.min(lastClickedIndex, currentIndex);
            const end = Math.max(lastClickedIndex, currentIndex);
            const rangeIndices: number[] = [];
            for (let i = 0; i < filteredIndices.length; i++) {
                const idx = filteredIndices[i];
                if (idx >= start && idx <= end) {
                    rangeIndices.push(idx);
                }
            }
            highlightRows(rangeIndices);
            setLastClickedIndex(currentIndex);
            return;
        }

        // Ctrl/Cmd-click: toggle this gene's rows in the current selection
        if (isToggle) {
            const currentSet = new Set(highlightedIndices);
            const allSelected = geneRowIndices.every((idx) => currentSet.has(idx));

            if (allSelected) {
                geneRowIndices.forEach((idx) => currentSet.delete(idx));
            } else {
                geneRowIndices.forEach((idx) => currentSet.add(idx));
            }

            const newIndices = Array.from(currentSet).sort((a, b) => a - b);
            highlightRows(newIndices);
            setLastClickedIndex(Math.min(...geneRowIndices));
            return;
        }

        // Plain click: replace selection with this gene's rows
        highlightRows(geneRowIndices);
        setLastClickedIndex(Math.min(...geneRowIndices));
    };

    return (
        <div className="absolute w-[100%] h-[100%] overflow-y-auto text-sm h-full flex flex-col">
            {visibleGeneIds.length > 0 ? (
                <>
                    {hasMore && (
                        <div className="p-3 mb-2 text-xs text-gray-500">
                            Showing first {maxGenes} of {geneIds.length} genes (adjust cap in settings)
                        </div>
                    )}
                    <div ref={scrollContainerRef} className="flex-1 min-h-0 overflow-y-auto">
                        {visibleGeneIds.map((geneId) => {
                            const highlightCount = geneHighlightCounts.get(geneId) || 0;
                            const isHighlighted = highlightCount > 0;
                            return (
                                <div
                                    key={geneId}
                                    ref={(el) => setGeneRef(geneId, el)}
                                    data-gene-id={geneId}
                                >
                                    <div className="relative">
                                        <GeneNetworkInfoComponent
                                            geneId={geneId}
                                            highlightCount={highlightCount}
                                            isHighlighted={isHighlighted}
                                            onCardClick={(event) => handleGeneClick(geneId, event)}
                                        />
                                    </div>
                                </div>
                            );
                        })}
                    </div>
                </>
            ) : (
                <div>
                    {mode === "observableFields" && highlightedIndices.length === 0
                        ? "No filtered data. Apply filters to see genes."
                        : "No gene IDs found."}
                </div>
            )}
        </div>
    );
});

