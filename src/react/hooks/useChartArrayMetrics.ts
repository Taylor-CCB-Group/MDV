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
        const root = layoutRef.current?.root;
        if (!root || cellCount === 0) {
            setMetrics(EMPTY_METRICS);
            return;
        }

        const measure = () => {
            const cells = layoutRef.current?.cells ?? [];
            setMetrics({
                rootSize: measureRootSize(root),
                cellBounds: measureCellRectsRelativeToRoot(root, cells.slice(0, cellCount)),
            });
        };

        measure();
        const observer = new ResizeObserver(measure);
        observer.observe(root);
        const cells = layoutRef.current?.cells ?? [];
        for (let i = 0; i < cellCount; i++) {
            const cell = cells[i];
            if (cell) observer.observe(cell);
        }
        window.addEventListener("resize", measure);
        return () => {
            observer.disconnect();
            window.removeEventListener("resize", measure);
        };
    }, [layoutRef, cellCount]);

    return metrics;
}
