import { useMemo, useState, useRef, useCallback } from "react";
import { observer } from "mobx-react-lite";
import { useChart } from "../context";
import { useParamColumns } from "../hooks";
import { useHighlightedIndex } from "../selectionHooks";
import { useSimplerFilteredIndices } from "../hooks";
import { isColumnLoaded } from "@/lib/columnTypeHelpers";
import { GeneNetworkInfoComponent } from "./GeneNetworkInfoComponent";
import type { GeneNetworkConfig } from "./GeneNetworkChart";

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

/**
 * Simple virtualized list component for displaying filtered gene info.
 * Only renders items that are visible in the viewport.
 */
function VirtualizedGeneList({ 
    geneIds, 
    itemHeight = 200 
}: { 
    geneIds: string[];
    itemHeight?: number;
}) {
    const [scrollTop, setScrollTop] = useState(0);
    const containerRef = useRef<HTMLDivElement>(null);

    const handleScroll = useCallback((e: React.UIEvent<HTMLDivElement>) => {
        setScrollTop(e.currentTarget.scrollTop);
    }, []);

    const visibleRange = useMemo(() => {
        const containerHeight = containerRef.current?.clientHeight || 600;
        const start = Math.floor(scrollTop / itemHeight);
        const end = Math.min(
            geneIds.length,
            Math.ceil((scrollTop + containerHeight) / itemHeight) + 1,
        );
        return { start, end };
    }, [scrollTop, itemHeight, geneIds.length]);

    const visibleItems = useMemo(() => {
        return geneIds.slice(visibleRange.start, visibleRange.end).map((geneId, i) => ({
            geneId,
            index: visibleRange.start + i,
        }));
    }, [geneIds, visibleRange]);

    const totalHeight = geneIds.length * itemHeight;
    const offsetY = visibleRange.start * itemHeight;

    return (
        <div
            ref={containerRef}
            className="h-full overflow-y-auto"
            onScroll={handleScroll}
        >
            <div style={{ height: `${totalHeight}px`, position: "relative" }}>
                <div
                    style={{
                        transform: `translateY(${offsetY}px)`,
                        position: "absolute",
                        width: "100%",
                    }}
                >
                    {visibleItems.map(({ geneId, index }) => (
                        <div
                            key={`${geneId}-${index}`}
                            style={{
                                height: `${itemHeight}px`,
                                padding: "0.5em",
                            }}
                        >
                            <GeneNetworkInfoComponent geneId={geneId} />
                        </div>
                    ))}
                </div>
            </div>
        </div>
    );
}

export const GeneNetworkChartComponent = observer(() => {
    const chart = useChart<GeneNetworkConfig>();
    const columns = useParamColumns();
    const highlightedIndex = useHighlightedIndex();
    const filteredIndices = useSimplerFilteredIndices();

    if (columns.length === 0) {
        return <div>No column selected. Please configure the chart.</div>;
    }

    const column = columns[0];
    if (!isColumnLoaded(column)) {
        return <div>Loading column data...</div>;
    }

    const displayMode = chart.config.displayMode || "highlighted";

    // Get gene IDs based on display mode
    const geneIds = useMemo(() => {
        if (displayMode === "highlighted") {
            if (highlightedIndex >= 0) {
                try {
                    const geneId = getTextValue(column, highlightedIndex);
                    return geneId ? [geneId] : [];
                } catch (e) {
                    console.error("Error getting highlighted gene ID:", e);
                    return [];
                }
            }
            return [];
        } else {
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
            return Array.from(geneIdSet);
        }
    }, [displayMode, highlightedIndex, filteredIndices, column]);

    return (
        <div className="absolute w-[100%] h-[100%] overflow-y-auto text-sm h-full flex flex-col">
            {displayMode === "highlighted" ? (
                <>
                    {highlightedIndex >= 0 ? (
                        <>
                            {geneIds.length > 0 ? (
                                <GeneNetworkInfoComponent geneId={geneIds[0]} />
                            ) : (
                                <div>No gene ID found for highlighted row.</div>
                            )}
                        </>
                    ) : (
                        <div>No row highlighted. Click on a data point to highlight it.</div>
                    )}
                </>
            ) : (
                <>
                    {geneIds.length > 0 ? (
                        <div className="flex-1 min-h-0">
                            <VirtualizedGeneList geneIds={geneIds} />
                        </div>
                    ) : (
                        <div>No gene IDs found in filtered data.</div>
                    )}
                </>
            )}
        </div>
    );
});

