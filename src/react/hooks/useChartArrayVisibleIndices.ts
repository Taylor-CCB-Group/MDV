import { useEffect, useState, type RefObject } from "react";
import type { ChartArrayLayoutHandle } from "../components/ChartArrayLayout";

export function useChartArrayVisibleIndices(
    scrollRootRef: RefObject<HTMLElement | null>,
    layoutRef: RefObject<ChartArrayLayoutHandle | null>,
    cellCount: number,
): number[] {
    const [visibleIndices, setVisibleIndices] = useState<number[]>(() =>
        cellCount > 0 ? Array.from({ length: cellCount }, (_, index) => index) : [],
    );

    useEffect(() => {
        const scrollRoot = scrollRootRef.current;
        const cells = layoutRef.current?.cells ?? [];
        if (!scrollRoot || cellCount === 0) {
            setVisibleIndices([]);
            return;
        }

        const intersectionByIndex = new Map<number, boolean>();

        const syncVisible = () => {
            const indices = [...intersectionByIndex.entries()]
                .filter(([, isVisible]) => isVisible)
                .map(([index]) => index)
                .sort((a, b) => a - b);
            setVisibleIndices(indices.length > 0 ? indices : []);
        };

        const observer = new IntersectionObserver(
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

        for (let i = 0; i < cellCount; i++) {
            const cell = cells[i];
            if (cell) observer.observe(cell);
        }

        return () => observer.disconnect();
    }, [scrollRootRef, layoutRef, cellCount]);

    return visibleIndices;
}
