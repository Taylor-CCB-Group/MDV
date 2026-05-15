import { useLayoutEffect, useState, type RefObject } from "react";
import type { ChartArrayLayoutHandle } from "../components/ChartArrayLayout";

function allCellIndices(cellCount: number) {
    return Array.from({ length: cellCount }, (_, index) => index);
}

export function useChartArrayVisibleIndices(
    scrollRootRef: RefObject<HTMLElement | null>,
    layoutRef: RefObject<ChartArrayLayoutHandle | null>,
    cellCount: number,
): number[] {
    const [visibleIndices, setVisibleIndices] = useState<number[]>(() =>
        cellCount > 0 ? allCellIndices(cellCount) : [],
    );

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
                if (setupAttempts < 20) {
                    setupAttempts += 1;
                    requestAnimationFrame(setup);
                }
                return;
            }

            if (!observer) {
                observer = new IntersectionObserver(
                    (entries) => {
                        for (const entry of entries) {
                            const index = Number((entry.target as HTMLElement).dataset.chartArrayIndex);
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
    }, [scrollRootRef, layoutRef, cellCount]);

    return visibleIndices;
}
