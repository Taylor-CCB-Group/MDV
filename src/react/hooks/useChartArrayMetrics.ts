import { type RefObject, useLayoutEffect, useState } from "react";
import type { ChartArrayLayoutHandle } from "../components/ChartArrayLayout";
import { type ChartArrayRect, measureCellRectsRelativeToRoot, measureRootSize } from "../components/chartArrayUtils";

export type ChartArrayMetrics = {
    rootSize: ChartArrayRect;
    cellBounds: ChartArrayRect[];
};

const EMPTY_METRICS: ChartArrayMetrics = {
    rootSize: { x: 0, y: 0, width: 0, height: 0 },
    cellBounds: [],
};

/**
 * The size of the layout and of each cell, kept current by a ResizeObserver. `cellsKey` names the
 * cells (their keys joined): a new cell in a slot is a new element, so a change of key observes
 * the cells again, even when their number is the same.
 */
export function useChartArrayMetrics(
    layoutRef: RefObject<ChartArrayLayoutHandle | null>,
    cellCount: number,
    cellsKey: string,
): ChartArrayMetrics {
    const [metrics, setMetrics] = useState<ChartArrayMetrics>(EMPTY_METRICS);

    // biome-ignore lint/correctness/useExhaustiveDependencies: cellsKey is the trigger for new cell elements
    useLayoutEffect(() => {
        if (cellCount === 0) {
            setMetrics(EMPTY_METRICS);
            return;
        }

        let cancelled = false;
        let resizeObserver: ResizeObserver | null = null;

        const measure = () => {
            if (cancelled) return;
            const root = layoutRef.current?.root;
            if (!root) return;
            const cells = layoutRef.current?.cells ?? [];
            setMetrics({
                rootSize: measureRootSize(root),
                cellBounds: measureCellRectsRelativeToRoot(root, cells.slice(0, cellCount)),
            });
        };

        const setup = () => {
            if (cancelled) return;
            const root = layoutRef.current?.root;
            if (!root) {
                // the layout can mount well after this (a grid shows a placeholder until its
                // device is ready), and nothing else would set this up then
                requestAnimationFrame(setup);
                return;
            }

            resizeObserver?.disconnect();
            measure();
            resizeObserver = new ResizeObserver(measure);
            resizeObserver.observe(root);
            const cells = layoutRef.current?.cells ?? [];
            for (let i = 0; i < cellCount; i++) {
                const cell = cells[i];
                if (cell) resizeObserver.observe(cell);
            }
        };

        setup();
        window.addEventListener("resize", measure);

        return () => {
            cancelled = true;
            resizeObserver?.disconnect();
            window.removeEventListener("resize", measure);
        };
    }, [layoutRef, cellCount, cellsKey]);

    return metrics;
}
