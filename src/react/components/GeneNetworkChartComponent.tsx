import { useMemo, useRef, useLayoutEffect, useState, useEffect, useCallback, type MouseEvent, type KeyboardEvent } from "react";
import { observer } from "mobx-react-lite";
import { useVirtualizer } from "@tanstack/react-virtual";
import { useChart } from "../context";
import { useParamColumns } from "../hooks";
import { useHighlightedIndices, useHighlightRows } from "../selectionHooks";
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

export const GeneNetworkChartComponent = observer(() => {
    const chart = useChart<GeneNetworkConfig>();
    const [column] = useParamColumns();
    const highlightedIndices = useHighlightedIndices();
    const filteredIndices = useSimplerFilteredIndices();
    const highlightRows = useHighlightRows();
    const scrollContainerRef = useRef<HTMLDivElement>(null);
    const geneRefs = useRef<Map<string, HTMLDivElement>>(new Map());
    const [lastClickedIndex, setLastClickedIndex] = useState<number | null>(null);
    const [focusedGeneIndex, setFocusedGeneIndex] = useState<number | null>(null);
    const [rangeAnchorGeneIndex, setRangeAnchorGeneIndex] = useState<number | null>(null);
    const prevHighlightedIndicesRef = useRef<number[] | null>(null);

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

    // Set up virtualizer for efficient rendering of large gene lists
    // Using fixed height of 180px (matches GeneNetworkInfoComponent fixed height)
    const rowVirtualizer = useVirtualizer({
        count: visibleGeneIds.length,
        getScrollElement: () => scrollContainerRef.current,
        estimateSize: () => 180, // Fixed height of gene card (matches GeneNetworkInfoComponent)
        overscan: 5, // Render a few extra items for smoother scrolling
    });

    // Auto-scroll to highlighted genes from external sources (smart interaction detection)
    // Default to true unless explicitly disabled in config
    const autoScroll = chart.config.autoScroll ?? true;
    useLayoutEffect(() => {
        if (!autoScroll) return;
        if (highlightedIndices.length === 0 || visibleGeneIds.length === 0) return;
        if (!scrollContainerRef.current) return;

        // Only auto-scroll when highlighted indices actually change (not on scroll/layout changes)
        const prevHighlighted = prevHighlightedIndicesRef.current;
        const highlightedChanged = !prevHighlighted || 
            prevHighlighted.length !== highlightedIndices.length ||
            !prevHighlighted.every((val, idx) => val === highlightedIndices[idx]);
        
        if (!highlightedChanged) {
            return; // Highlighted indices haven't changed, don't auto-scroll
        }
        
        // Update ref for next comparison
        prevHighlightedIndicesRef.current = highlightedIndices;

        // Check if user is actively interacting with the component
        // If so, don't auto-scroll to avoid feedback loops
        const hasComponentFocus = scrollContainerRef.current?.contains(document.activeElement) ?? false;
        const isUserInteracting = hasComponentFocus && (
            focusedGeneIndex !== null || // User navigating with keyboard
            rangeAnchorGeneIndex !== null // User extending range selection
        );

        if (isUserInteracting) {
            return; // Don't auto-scroll while user is interacting with this component
        }

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

        // Find the first highlighted gene that's outside the actual viewport
        // Calculate visible range based on scroll position and container height (not overscan)
        const container = scrollContainerRef.current;
        if (!container) return;
        
        const scrollTop = container.scrollTop;
        const containerHeight = container.clientHeight;
        const itemHeight = 180; // Fixed height from estimateSize
        
        // Calculate which indices are actually visible in the viewport
        const firstVisibleIndex = Math.floor(scrollTop / itemHeight);
        const lastVisibleIndex = Math.min(
            visibleGeneIds.length - 1,
            Math.floor((scrollTop + containerHeight) / itemHeight)
        );
        
        // Find first highlighted gene that's outside the actual viewport
        for (let i = 0; i < visibleGeneIds.length; i++) {
            const geneId = visibleGeneIds[i];
            if (highlightedGeneIds.has(geneId)) {
                // Check if this gene is outside the actual viewport (with small buffer)
                if (i < firstVisibleIndex || i > lastVisibleIndex) {
                    // Found a highlighted gene outside viewport - scroll to it
                    rowVirtualizer.scrollToIndex(i, {
                        align: "center",
                        behavior: "smooth",
                    });
                    break; // Only scroll to first one
                }
            }
        }
    }, [highlightedIndices, visibleGeneIds, column, autoScroll, focusedGeneIndex, rangeAnchorGeneIndex, rowVirtualizer]);

    const setGeneRef = (geneId: string, element: HTMLDivElement | null) => {
        if (element) {
            geneRefs.current.set(geneId, element);
        } else {
            geneRefs.current.delete(geneId);
        }
    };

    const getIndicesForGeneId = useCallback((geneId: string): number[] => {
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
    }, [filteredIndices, column]);

    const handleGeneSelection = useCallback((
        geneId: string,
        isToggle: boolean,
        isRange: boolean,
        currentGeneIndex?: number,
    ) => {
        // including a reference so linter doesn't think this is an unneeded dependency
        filteredIndices; // changes to this are important, even though not used directly
        const geneRowIndices = getIndicesForGeneId(geneId);
        if (geneRowIndices.length === 0) return;

        // Shift: select a continuous range of genes (not just row indices)
        if (isRange && currentGeneIndex !== undefined) {
            // Use rangeAnchorGeneIndex if set, otherwise establish a new anchor
            let anchorGeneIndex: number;
            if (rangeAnchorGeneIndex !== null) {
                anchorGeneIndex = rangeAnchorGeneIndex;
            } else {
                // Starting a new range selection - find anchor from lastClickedIndex or use current
                anchorGeneIndex = -1;
                if (lastClickedIndex !== null) {
                    // Find which gene corresponds to lastClickedIndex
                    for (let i = 0; i < visibleGeneIds.length; i++) {
                        const testGeneId = visibleGeneIds[i];
                        const testIndices = getIndicesForGeneId(testGeneId);
                        if (testIndices.some(idx => idx === lastClickedIndex)) {
                            anchorGeneIndex = i;
                            break;
                        }
                    }
                }
                // If we couldn't find it, use the current gene index as anchor
                if (anchorGeneIndex === -1) {
                    anchorGeneIndex = currentGeneIndex;
                }
                setRangeAnchorGeneIndex(anchorGeneIndex);
            }

            // Select all genes from anchor to current
            const start = Math.min(anchorGeneIndex, currentGeneIndex);
            const end = Math.max(anchorGeneIndex, currentGeneIndex);
            const allRangeIndices: number[] = [];
            for (let i = start; i <= end; i++) {
                const rangeGeneId = visibleGeneIds[i];
                const rangeGeneIndices = getIndicesForGeneId(rangeGeneId);
                allRangeIndices.push(...rangeGeneIndices);
            }
            highlightRows(allRangeIndices);
            const currentIndex = Math.min(...geneRowIndices);
            setLastClickedIndex(currentIndex);
            return;
        }

        // Ctrl/Cmd: toggle this gene's rows in the current selection
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

        // Plain: replace selection with this gene's rows
        highlightRows(geneRowIndices);
        setLastClickedIndex(Math.min(...geneRowIndices));
        // Clear range anchor when doing plain selection
        setRangeAnchorGeneIndex(null);
    }, [filteredIndices, highlightedIndices, lastClickedIndex, highlightRows, getIndicesForGeneId, visibleGeneIds, rangeAnchorGeneIndex]);

    const focusGeneAtIndex = useCallback((index: number, shouldSelect: boolean, isRange: boolean, isToggle: boolean) => {
        if (index < 0 || index >= visibleGeneIds.length) return;
        
        const geneId = visibleGeneIds[index];
        setFocusedGeneIndex(index);
        
        // Scroll to the gene using virtualizer
        rowVirtualizer.scrollToIndex(index, {
            align: "center",
            behavior: "smooth",
        });
        
        // Focus the gene card element
        const geneElement = geneRefs.current.get(geneId);
        if (geneElement) {
            // Find the Paper element (GeneNetworkInfoComponent) within the gene element
            const paperElement = geneElement.querySelector('[tabindex="0"]') as HTMLElement;
            if (paperElement) {
                paperElement.focus();
            }
        }
        
        if (shouldSelect) {
            handleGeneSelection(geneId, isToggle, isRange, index);
        }
    }, [visibleGeneIds, rowVirtualizer, handleGeneSelection]);

    // Container-level keyboard handler for arrow key navigation
    const handleContainerKeyDown = useCallback((event: KeyboardEvent<HTMLDivElement>) => {
        // Only handle arrow keys
        if (event.key !== "ArrowUp" && event.key !== "ArrowDown") {
            return;
        }

        // Don't handle if user is typing in an input/textarea
        const target = event.target as HTMLElement;
        if (target.tagName === "INPUT" || target.tagName === "TEXTAREA" || target.isContentEditable) {
            return;
        }

        event.preventDefault();
        event.stopPropagation();

        const isRange = event.shiftKey;
        const isToggle = event.ctrlKey || event.metaKey;
        const shouldSelect = !isToggle; // Ctrl/Cmd+Arrow just navigates without selecting

        let newIndex: number;
        if (event.key === "ArrowDown") {
            if (focusedGeneIndex === null) {
                newIndex = 0;
            } else {
                newIndex = Math.min(focusedGeneIndex + 1, visibleGeneIds.length - 1);
            }
        } else {
            // ArrowUp
            if (focusedGeneIndex === null) {
                newIndex = visibleGeneIds.length - 1;
            } else {
                newIndex = Math.max(focusedGeneIndex - 1, 0);
            }
        }

        focusGeneAtIndex(newIndex, shouldSelect, isRange, isToggle);
    }, [focusedGeneIndex, visibleGeneIds, focusGeneAtIndex]);

    const handleGeneClick = (geneId: string, event: MouseEvent<HTMLDivElement>) => {
        event.preventDefault();
        event.stopPropagation();

        // Update focused gene index when clicked
        const index = visibleGeneIds.indexOf(geneId);
        if (index !== -1) {
            setFocusedGeneIndex(index);
        }

        const isToggle = event.ctrlKey || event.metaKey;
        const isRange = event.shiftKey;
        
        // Clear range anchor when clicking (not extending with keyboard)
        if (!isRange) {
            setRangeAnchorGeneIndex(null);
        }
        
        handleGeneSelection(geneId, isToggle, isRange, index);
    };

    const handleGeneKeyDown = (geneId: string, event: KeyboardEvent<HTMLDivElement>) => {
        // Handle arrow keys - let them bubble to container handler
        if (event.key === "ArrowUp" || event.key === "ArrowDown") {
            // Update focused index when arrow key is pressed on a gene card
            const index = visibleGeneIds.indexOf(geneId);
            if (index !== -1) {
                setFocusedGeneIndex(index);
            }
            // Let the event bubble to container handler for navigation
            return;
        }

        // Only handle Space key for selection
        if (event.key !== " " && event.key !== "Spacebar") {
            return;
        }

        // Avoid interfering with OS-level shortcuts (e.g. meta+space)
        if (event.metaKey || event.ctrlKey) {
            return;
        }

        event.preventDefault();
        event.stopPropagation();

        // Update focused gene index when Space is pressed
        const index = visibleGeneIds.indexOf(geneId);
        if (index !== -1) {
            setFocusedGeneIndex(index);
        }

        // Keyboard Space:
        // - plain Space toggles this gene
        // - Shift+Space selects a range (like Shift-click)
        const isRange = event.shiftKey;
        const isToggle = !isRange;
        
        // Clear range anchor when using Space (not extending with keyboard)
        if (!isRange) {
            setRangeAnchorGeneIndex(null);
        }
        
        handleGeneSelection(geneId, isToggle, isRange, index);
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
                    <div 
                        ref={scrollContainerRef} 
                        className="flex-1 min-h-0 overflow-y-auto"
                        onKeyDown={handleContainerKeyDown}
                    >
                        <div
                            style={{
                                height: `${rowVirtualizer.getTotalSize()}px`,
                                width: "100%",
                                position: "relative",
                            }}
                        >
                            {rowVirtualizer.getVirtualItems().map((virtualItem: { index: number; start: number }) => {
                                const geneId = visibleGeneIds[virtualItem.index];
                                const highlightCount = geneHighlightCounts.get(geneId) || 0;
                                const isHighlighted = highlightCount > 0;
                                return (
                                    <div
                                        key={geneId}
                                        ref={(el) => setGeneRef(geneId, el)}
                                        data-gene-id={geneId}
                                        data-index={virtualItem.index}
                                        style={{
                                            position: "absolute",
                                            top: `${virtualItem.start}px`,
                                            left: 0,
                                            width: "100%",
                                        }}
                                    >
                                        <div className="relative">
                                            <GeneNetworkInfoComponent
                                                geneId={geneId}
                                                isHighlighted={isHighlighted}
                                                onCardClick={(event) => handleGeneClick(geneId, event)}
                                                onCardKeyDown={(event) => handleGeneKeyDown(geneId, event)}
                                            />
                                        </div>
                                    </div>
                                );
                            })}
                        </div>
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

