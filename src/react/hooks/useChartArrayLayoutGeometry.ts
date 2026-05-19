/**
 * Resize-driven column count and viewport-fill detection for chart-array grid mode.
 * Prototype implementation — likely to be replaced or simplified when layout becomes user-configurable.
 */
import { useLayoutEffect, useState, type RefObject } from "react";
import {
    chartArrayGridFitsViewport,
    getPreferredGridColumnCount,
    readChartArraySizingFromRoot,
    type ChartArrayLayoutMode,
} from "../components/chartArrayLayoutUtils";

export type ChartArrayLayoutGeometry = {
    gridFillsViewport: boolean;
    gridColumns: number;
};

const DEFAULT_GEOMETRY: ChartArrayLayoutGeometry = {
    gridFillsViewport: true,
    gridColumns: 1,
};

export function useChartArrayLayoutGeometry(
    layoutRootRef: RefObject<HTMLDivElement | null>,
    cellCount: number,
    layoutMode: ChartArrayLayoutMode,
): ChartArrayLayoutGeometry {
    const [geometry, setGeometry] = useState(DEFAULT_GEOMETRY);

    useLayoutEffect(() => {
        if (layoutMode !== "grid" || cellCount === 0) {
            setGeometry(DEFAULT_GEOMETRY);
            return;
        }

        let cancelled = false;
        let resizeObserver: ResizeObserver | null = null;

        const measure = () => {
            if (cancelled) return;
            const root = layoutRootRef.current;
            const scrollRoot = root?.parentElement;
            if (!root || !scrollRoot) return;

            const sizing = readChartArraySizingFromRoot(root);
            const innerWidth = Math.max(0, root.clientWidth - sizing.paddingX);
            const viewportHeight = scrollRoot.clientHeight;
            setGeometry({
                gridColumns: getPreferredGridColumnCount(cellCount, innerWidth, sizing),
                gridFillsViewport: chartArrayGridFitsViewport(
                    cellCount,
                    viewportHeight,
                    innerWidth,
                    sizing,
                ),
            });
        };

        let setupAttempts = 0;
        const setup = () => {
            if (cancelled) return;
            const root = layoutRootRef.current;
            const scrollRoot = root?.parentElement;
            if (!root || !scrollRoot) {
                if (setupAttempts < 20) {
                    setupAttempts += 1;
                    requestAnimationFrame(setup);
                }
                return;
            }

            resizeObserver?.disconnect();
            measure();
            resizeObserver = new ResizeObserver(measure);
            resizeObserver.observe(scrollRoot);
            resizeObserver.observe(root);
        };

        setup();
        window.addEventListener("resize", measure);

        return () => {
            cancelled = true;
            resizeObserver?.disconnect();
            window.removeEventListener("resize", measure);
        };
    }, [layoutRootRef, cellCount, layoutMode]);

    if (layoutMode !== "grid") {
        return DEFAULT_GEOMETRY;
    }
    return geometry;
}
