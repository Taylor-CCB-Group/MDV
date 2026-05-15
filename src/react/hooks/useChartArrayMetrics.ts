import { useLayoutEffect, useState, type RefObject } from "react";
import type { ChartArrayLayoutHandle } from "../components/ChartArrayLayout";
import {
    measureCellRectsRelativeToRoot,
    measureRootSize,
    type ChartArrayRect,
} from "../components/chartArrayUtils";

export type ChartArrayMetrics = {
    rootSize: ChartArrayRect;
    cellBounds: ChartArrayRect[];
};

const EMPTY_METRICS: ChartArrayMetrics = {
    rootSize: { x: 0, y: 0, width: 0, height: 0 },
    cellBounds: [],
};

export function useChartArrayMetrics(
    layoutRef: RefObject<ChartArrayLayoutHandle | null>,
    cellCount: number,
): ChartArrayMetrics {
    const [metrics, setMetrics] = useState<ChartArrayMetrics>(EMPTY_METRICS);

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

        let setupAttempts = 0;
        const setup = () => {
            if (cancelled) return;
            const root = layoutRef.current?.root;
            if (!root) {
                if (setupAttempts < 20) {
                    setupAttempts += 1;
                    requestAnimationFrame(setup);
                }
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
    }, [layoutRef, cellCount]);

    return metrics;
}
