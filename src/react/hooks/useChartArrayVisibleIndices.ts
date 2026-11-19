import { type RefObject, useLayoutEffect, useState } from "react";
import type { ChartArrayLayoutHandle } from "../components/ChartArrayLayout";

function allCellIndices(cellCount: number) {
    return Array.from({ length: cellCount }, (_, index) => index);
}

/**
 * The indices of the cells in or near view, from an IntersectionObserver. `cellsKey` names the
 * cells (their keys joined): a new cell in a slot is a new element, so a change of key observes
 * the cells again, from a clean slate, even when their number is the same.
 */
export function useChartArrayVisibleIndices(
    scrollRootRef: RefObject<HTMLElement | null>,
    layoutRef: RefObject<ChartArrayLayoutHandle | null>,
    cellCount: number,
    cellsKey: string,
): number[] {
    const [visibleIndices, setVisibleIndices] = useState<number[]>(() =>
        cellCount > 0 ? allCellIndices(cellCount) : [],
    );

    // biome-ignore lint/correctness/useExhaustiveDependencies: cellsKey is the trigger for new cell elements
    useLayoutEffect(() => {
        if (cellCount === 0) {
            setVisibleIndices([]);
            return;
        }

        let cancelled = false;
        let observer: IntersectionObserver | null = null;
        let setupAttempts = 0;
        const intersectionByIndex = new Map<number, boolean>();
        const observedCells = new WeakSet<HTMLElement>();

        const syncVisible = () => {
            const indices = [...intersectionByIndex.entries()]
                .filter(([, isVisible]) => isVisible)
                .map(([index]) => index)
                .sort((a, b) => a - b);
            setVisibleIndices(indices.length > 0 ? indices : allCellIndices(cellCount));
        };

        const setup = () => {
            if (cancelled) return;
            const scrollRoot = scrollRootRef.current;
            if (!scrollRoot) {
                // the grid can mount well after this (it shows a placeholder until its device
                // is ready), and nothing else would set this up then
                requestAnimationFrame(setup);
                return;
            }

            if (!observer) {
                observer = new IntersectionObserver(
                    (entries) => {
                        for (const entry of entries) {
                            // not `instanceof HTMLElement`: a cell made in a popout's document
                            // is that window's, and would fail it
                            const index = Number(entry.target.getAttribute("data-chart-array-index"));
                            if (!Number.isInteger(index) || index < 0 || index >= cellCount) continue;
                            intersectionByIndex.set(index, entry.isIntersecting);
                        }
                        syncVisible();
                    },
                    { root: scrollRoot, rootMargin: "64px" },
                );
            }

            const cells = layoutRef.current?.cells ?? [];
            let observed = 0;
            for (let i = 0; i < cellCount; i++) {
                const cell = cells[i];
                if (cell) {
                    if (!observedCells.has(cell)) {
                        observer.observe(cell);
                        observedCells.add(cell);
                    }
                    observed += 1;
                }
            }
            syncVisible();
            if (observed < cellCount && setupAttempts < 20) {
                setupAttempts += 1;
                requestAnimationFrame(setup);
            }
        };

        setup();

        return () => {
            cancelled = true;
            observer?.disconnect();
        };
    }, [scrollRootRef, layoutRef, cellCount, cellsKey]);

    return visibleIndices;
}
