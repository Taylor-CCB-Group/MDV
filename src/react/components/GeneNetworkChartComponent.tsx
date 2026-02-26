import { useMemo, useRef, useLayoutEffect, useState, useCallback, type MouseEvent, type KeyboardEvent } from "react";
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
    const rowRefs = useRef<Map<number, HTMLDivElement>>(new Map());
    const [lastClickedIndex, setLastClickedIndex] = useState<number | null>(null);
    const [focusedGeneIndex, setFocusedGeneIndex] = useState<number | null>(null);
    const [rangeAnchorGeneIndex, setRangeAnchorGeneIndex] = useState<number | null>(null);
    const prevHighlightedIndicesRef = useRef<number[] | null>(null);
    const lastScrolledGeneIndexRef = useRef<number | null>(null);

    // Row-centric model: visible list is filtered indices only; no string identity, no getTextValue in memo.
    const { visibleRowIndices, rowIndexToVisibleIndex, highlightedRowSet } = useMemo(() => {
        const visible = [...filteredIndices];
        const map = new Map<number, number>();
        for (let i = 0; i < visible.length; i++) {
            map.set(visible[i], i);
        }
        const highlighted = new Set(highlightedIndices);
        return {
            visibleRowIndices: visible,
            rowIndexToVisibleIndex: map,
            highlightedRowSet: highlighted,
        };
    }, [filteredIndices, highlightedIndices]);

    const ROW_HEIGHT = 150;
    const ROW_VERTICAL_GAP = 10; // vertical space between items
    const ROW_TOTAL_HEIGHT = ROW_HEIGHT + ROW_VERTICAL_GAP;

    // Set up virtualizer for efficient rendering of large gene lists
    // Using fixed height of gene card plus a vertical gap between items
    const rowVirtualizer = useVirtualizer({
        count: visibleRowIndices.length,
        getScrollElement: () => scrollContainerRef.current,
        estimateSize: () => ROW_TOTAL_HEIGHT,
        overscan: 5, // Render a few extra items for smoother scrolling
    });

    // Auto-scroll to highlighted genes from external sources (smart interaction detection).
    // Controlled via chart config so it can be toggled from settings and serialised.
    const autoScroll = chart.config.autoScroll ?? true;
    useLayoutEffect(() => {
        if (!autoScroll) return;
        if (highlightedIndices.length === 0 || visibleRowIndices.length === 0) return;
        if (!scrollContainerRef.current) return;

        // Only auto-scroll when highlighted indices actually change (not on scroll/layout changes)
        const prevHighlighted = prevHighlightedIndicesRef.current;
        const highlightedChanged =
            !prevHighlighted ||
            prevHighlighted.length !== highlightedIndices.length ||
            !prevHighlighted.every((val, idx) => val === highlightedIndices[idx]);

        if (!highlightedChanged) {
            return; // Highlighted indices haven't changed, don't auto-scroll
        }

        // Compute sorted previous and current highlights for diffing
        const prevSorted = prevHighlighted ? [...prevHighlighted].sort((a, b) => a - b) : [];
        const currSorted = [...highlightedIndices].sort((a, b) => a - b);

        const prevSet = new Set(prevSorted);
        const currSet = new Set(currSorted);

        const added: number[] = currSorted.filter((i) => !prevSet.has(i));
        const removed: number[] = prevSorted.filter((i) => !currSet.has(i));

        // Decide target row index based on how the selection changed
        let targetRowIndex: number | null = null;

        if (currSorted.length === 1) {
            // If new selection has a single item, always scroll to it
            targetRowIndex = currSorted[0];
        } else if (added.length > 0 && removed.length === 0) {
            // Pure extension (typical Shift+Arrow/Shift+Click) -> scroll to new far edge
            targetRowIndex = added[added.length - 1];
        } else if (removed.length > 0 && added.length === 0) {
            // Pure contraction: scroll to the corresponding end
            const prevMin = prevSorted[0];
            const prevMax = prevSorted[prevSorted.length - 1];
            const currMin = currSorted[0];
            const currMax = currSorted[currSorted.length - 1];

            // Determine which end contracted more by comparing distances
            const contractedFromTop = currMin > prevMin;
            const contractedFromBottom = currMax < prevMax;

            if (contractedFromBottom && !contractedFromTop) {
                // Contracted from bottom: scroll to new bottom
                targetRowIndex = currMax;
            } else if (contractedFromTop && !contractedFromBottom) {
                // Contracted from top: scroll to new top
                targetRowIndex = currMin;
            } else {
                // Mixed or unclear contraction: pick closer end to previous focus
                if (focusedGeneIndex !== null) {
                    const focusedRowIndex = visibleRowIndices[focusedGeneIndex];
                    let closest = currSorted[0];
                    let minDist = Math.abs(closest - focusedRowIndex);
                    for (const idx of currSorted) {
                        const dist = Math.abs(idx - focusedRowIndex);
                        if (dist < minDist) {
                            closest = idx;
                            minDist = dist;
                        }
                    }
                    targetRowIndex = closest;
                } else {
                    // Fallback: scroll to new top
                    targetRowIndex = currMin;
                }
            }
        } else if (added.length > 0) {
            // Mixed changes (external or complex update): scroll to first newly added
            targetRowIndex = added[0];
        } else {
            // No clear target; don't scroll
            prevHighlightedIndicesRef.current = highlightedIndices;
            return;
        }

        // Update ref for next comparison
        prevHighlightedIndicesRef.current = highlightedIndices;

        // Check if user is actively interacting with the component
        // If so, don't auto-scroll to avoid feedback loops
        const hasComponentFocus =
            scrollContainerRef.current?.contains(document.activeElement) ?? false;
        const isUserInteracting =
            hasComponentFocus &&
            (focusedGeneIndex !== null || rangeAnchorGeneIndex !== null);

        if (isUserInteracting) {
            return; // Don't auto-scroll while user is interacting with this component
        }

        const targetVisibleIndex = rowIndexToVisibleIndex.get(targetRowIndex) ?? null;
        if (targetVisibleIndex == null) {
            return;
        }

        // Calculate visible range based on scroll position and container height
        const container = scrollContainerRef.current;
        if (!container) return;

        const scrollTop = container.scrollTop;
        const containerHeight = container.clientHeight;

        const targetRowTop = targetVisibleIndex * ROW_TOTAL_HEIGHT;
        const targetRowBottom = targetRowTop + ROW_TOTAL_HEIGHT;
        const viewportBottom = scrollTop + containerHeight;
        const isFullyVisible =
            targetRowTop >= scrollTop && targetRowBottom <= viewportBottom;

        // Only scroll if the target row is not fully visible in the viewport
        if (!isFullyVisible) {
            if (lastScrolledGeneIndexRef.current !== targetVisibleIndex) {
                lastScrolledGeneIndexRef.current = targetVisibleIndex;
                rowVirtualizer.scrollToIndex(targetVisibleIndex, {
                    align: "center",
                    behavior: "auto",
                });
            }
        } else if (lastScrolledGeneIndexRef.current === targetVisibleIndex) {
            lastScrolledGeneIndexRef.current = null;
        }
    }, [
        highlightedIndices,
        visibleRowIndices,
        rowIndexToVisibleIndex,
        autoScroll,
        focusedGeneIndex,
        rangeAnchorGeneIndex,
        rowVirtualizer,
        ROW_TOTAL_HEIGHT
    ]);

    const setRowRef = (rowIndex: number, element: HTMLDivElement | null) => {
        if (element) {
            rowRefs.current.set(rowIndex, element);
        } else {
            rowRefs.current.delete(rowIndex);
        }
    };

    const handleGeneSelection = useCallback((
        rowIndex: number,
        isToggle: boolean,
        isRange: boolean,
        currentVisibleIndex?: number,
    ) => {
        if (isRange && currentVisibleIndex !== undefined) {
            const maxVisibleIndex = visibleRowIndices.length - 1;
            if (maxVisibleIndex < 0) return;

            let anchorVisibleIndex: number;
            if (rangeAnchorGeneIndex !== null) {
                anchorVisibleIndex = Math.min(
                    Math.max(rangeAnchorGeneIndex, 0),
                    maxVisibleIndex,
                );
            } else {
                if (lastClickedIndex !== null) {
                    const anchorFromMap = rowIndexToVisibleIndex.get(lastClickedIndex);
                    anchorVisibleIndex = anchorFromMap ?? currentVisibleIndex;
                } else {
                    anchorVisibleIndex = currentVisibleIndex;
                }
                setRangeAnchorGeneIndex(anchorVisibleIndex);
            }
            const safeCurrentVisibleIndex = Math.min(
                Math.max(currentVisibleIndex, 0),
                maxVisibleIndex,
            );
            const start = Math.min(anchorVisibleIndex, safeCurrentVisibleIndex);
            const end = Math.max(anchorVisibleIndex, safeCurrentVisibleIndex);
            const allRangeIndices: number[] = [];
            for (let i = start; i <= end; i++) {
                const nextRowIndex = visibleRowIndices[i];
                if (nextRowIndex !== undefined) allRangeIndices.push(nextRowIndex);
            }
            highlightRows(allRangeIndices);
            setLastClickedIndex(rowIndex);
            return;
        }

        if (isToggle) {
            const currentSet = new Set(highlightedIndices);
            if (currentSet.has(rowIndex)) {
                currentSet.delete(rowIndex);
            } else {
                currentSet.add(rowIndex);
            }
            const newIndices = Array.from(currentSet).sort((a, b) => a - b);
            highlightRows(newIndices);
            setLastClickedIndex(rowIndex);
            return;
        }

        highlightRows([rowIndex]);
        setLastClickedIndex(rowIndex);
        setRangeAnchorGeneIndex(null);
    }, [highlightedIndices, lastClickedIndex, highlightRows, visibleRowIndices, rowIndexToVisibleIndex, rangeAnchorGeneIndex]);

    const focusGeneAtIndex = useCallback((visibleIndex: number, shouldSelect: boolean, isRange: boolean, isToggle: boolean) => {
        if (visibleIndex < 0 || visibleIndex >= visibleRowIndices.length) return;
        const rowIndex = visibleRowIndices[visibleIndex];
        setFocusedGeneIndex(visibleIndex);

        rowVirtualizer.scrollToIndex(visibleIndex, {
            align: "center",
            behavior: "auto",
        });

        const rowElement = rowRefs.current.get(rowIndex);
        if (rowElement) {
            const paperElement = rowElement.querySelector('[tabindex="0"]') as HTMLElement;
            if (paperElement) paperElement.focus();
        }

        if (shouldSelect) {
            handleGeneSelection(rowIndex, isToggle, isRange, visibleIndex);
        }
    }, [visibleRowIndices, rowVirtualizer, handleGeneSelection]);

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
                newIndex = Math.min(focusedGeneIndex + 1, visibleRowIndices.length - 1);
            }
        } else {
            if (focusedGeneIndex === null) {
                newIndex = visibleRowIndices.length - 1;
            } else {
                newIndex = Math.max(focusedGeneIndex - 1, 0);
            }
        }

        focusGeneAtIndex(newIndex, shouldSelect, isRange, isToggle);
    }, [focusedGeneIndex, visibleRowIndices.length, focusGeneAtIndex]);

    const handleGeneClick = (rowIndex: number, event: MouseEvent<HTMLDivElement>) => {
        event.preventDefault();
        event.stopPropagation();
        const visibleIndex = rowIndexToVisibleIndex.get(rowIndex);
        if (visibleIndex !== undefined) setFocusedGeneIndex(visibleIndex);
        const isToggle = event.ctrlKey || event.metaKey;
        const isRange = event.shiftKey;
        if (!isRange) setRangeAnchorGeneIndex(null);
        handleGeneSelection(rowIndex, isToggle, isRange, visibleIndex);
    };

    const handleGeneKeyDown = (rowIndex: number, event: KeyboardEvent<HTMLDivElement>) => {
        if (event.key === "ArrowUp" || event.key === "ArrowDown") {
            const visibleIndex = rowIndexToVisibleIndex.get(rowIndex);
            if (visibleIndex !== undefined) setFocusedGeneIndex(visibleIndex);
            return;
        }
        if (event.key !== " " && event.key !== "Spacebar") return;
        if (event.metaKey || event.ctrlKey) return;
        event.preventDefault();
        event.stopPropagation();
        const visibleIndex = rowIndexToVisibleIndex.get(rowIndex);
        if (visibleIndex !== undefined) setFocusedGeneIndex(visibleIndex);
        const isRange = event.shiftKey;
        const isToggle = !isRange;
        if (!isRange) setRangeAnchorGeneIndex(null);
        handleGeneSelection(rowIndex, isToggle, isRange, visibleIndex);
    };

    return (
        <div className="absolute w-[100%] h-[100%] text-sm h-full flex flex-col">
            {visibleRowIndices.length > 0 ? (
                <>
                    <div
                        ref={scrollContainerRef}
                        className="flex-1 min-h-0 overflow-y-auto"
                        onKeyDown={handleContainerKeyDown}
                    >
                        <div
                            className="relative w-full"
                            style={{ height: `${rowVirtualizer.getTotalSize()}px` }}
                        >
                            {rowVirtualizer.getVirtualItems().map((virtualItem: { index: number; start: number }) => {
                                const rowIndex = visibleRowIndices[virtualItem.index];
                                let geneId: string;
                                try {
                                    geneId = getTextValue(column, rowIndex);
                                } catch {
                                    geneId = "";
                                }
                                const isHighlighted = highlightedRowSet.has(rowIndex);
                                return (
                                    <div
                                        key={rowIndex}
                                        ref={(el) => setRowRef(rowIndex, el)}
                                        data-row-index={rowIndex}
                                        data-index={virtualItem.index}
                                        className="absolute left-0 w-full px-2"
                                        style={{
                                            top: `${virtualItem.start}px`,
                                            height: `${ROW_TOTAL_HEIGHT}px`,
                                        }}
                                    >
                                        <div className="relative h-full flex items-center">
                                            <GeneNetworkInfoComponent
                                                geneId={geneId}
                                                isHighlighted={isHighlighted}
                                                onCardClick={(event) => handleGeneClick(rowIndex, event)}
                                                onCardKeyDown={(event) => handleGeneKeyDown(rowIndex, event)}
                                            />
                                        </div>
                                    </div>
                                );
                            })}
                        </div>
                    </div>
                </>
            ) : (
                <div>No gene IDs found.</div>
            )}
        </div>
    );
});

