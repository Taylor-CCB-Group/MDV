/**
 * Resize-driven column count and viewport-fill detection for chart-array grid mode.
 * Prototype implementation — likely to be replaced or simplified when layout becomes user-configurable.
 */
import { type RefObject, useLayoutEffect, useState } from "react";
import {
    type ChartArrayCellSize,
    type ChartArrayLayoutMode,
    chartArrayCellSize,
    chooseChartArrayGrid,
    readChartArraySizingFromRoot,
} from "../components/chartArrayLayoutUtils";

export type ChartArrayLayoutGeometry = {
    gridFillsViewport: boolean;
    gridColumns: number;
    /** Every grid cell's size, in whole device pixels (see `chartArrayCellSize`); null until measured. */
    cellSize: ChartArrayCellSize | null;
};

const DEFAULT_GEOMETRY: ChartArrayLayoutGeometry = {
    gridFillsViewport: true,
    gridColumns: 1,
    cellSize: null,
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
            // client sizes are rounded, possibly up: cells sized from them could overflow and wrap
            const width = Math.min(root.clientWidth, root.getBoundingClientRect().width);
            const height = Math.min(scrollRoot.clientHeight, scrollRoot.getBoundingClientRect().height);
            const innerWidth = Math.max(0, width - sizing.paddingX);
            const innerHeight = Math.max(0, height - sizing.paddingY);
            const choice = chooseChartArrayGrid(cellCount, innerWidth, innerHeight, sizing);
            // the chart's own window: a popped-out chart may be on a screen of another density
            const pixelRatio = root.ownerDocument.defaultView?.devicePixelRatio ?? 1;
            const cellSize = chartArrayCellSize(cellCount, choice, innerWidth, innerHeight, pixelRatio, sizing);
            setGeometry((previous) =>
                previous.gridColumns === choice.columns &&
                previous.gridFillsViewport === choice.fillsViewport &&
                previous.cellSize?.width === cellSize.width &&
                previous.cellSize?.height === cellSize.height
                    ? previous
                    : { gridColumns: choice.columns, gridFillsViewport: choice.fillsViewport, cellSize },
            );
        };

        // moving to a screen of another pixel density resizes nothing, but changes the rounding
        let pixelRatioQuery: MediaQueryList | null = null;
        const watchPixelRatio = () => {
            pixelRatioQuery?.removeEventListener("change", onPixelRatioChange);
            const view = layoutRootRef.current?.ownerDocument.defaultView;
            if (!view) return;
            pixelRatioQuery = view.matchMedia(`(resolution: ${view.devicePixelRatio}dppx)`);
            pixelRatioQuery.addEventListener("change", onPixelRatioChange);
        };
        function onPixelRatioChange() {
            measure();
            watchPixelRatio();
        }

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
            watchPixelRatio();
            resizeObserver = new ResizeObserver(measure);
            resizeObserver.observe(scrollRoot);
            resizeObserver.observe(root);
        };

        setup();
        window.addEventListener("resize", measure);

        return () => {
            cancelled = true;
            resizeObserver?.disconnect();
            pixelRatioQuery?.removeEventListener("change", onPixelRatioChange);
            window.removeEventListener("resize", measure);
        };
    }, [layoutRootRef, cellCount, layoutMode]);

    if (layoutMode !== "grid") {
        return DEFAULT_GEOMETRY;
    }
    return geometry;
}
