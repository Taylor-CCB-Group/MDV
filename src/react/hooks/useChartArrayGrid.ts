import { useMemo, useRef, type RefObject } from "react";
import type { ChartArrayLayoutHandle } from "../components/ChartArrayLayout";
import type { ChartArrayRect } from "../components/chartArrayUtils";
import { getChartArrayViewId } from "../components/chartArrayGridUtils";
import { useChartArrayMetrics, type ChartArrayMetrics } from "./useChartArrayMetrics";
import { useChartArrayVisibleIndices } from "./useChartArrayVisibleIndices";

export type ChartArrayCell = {
    key: string;
    label: string;
};

export type UseChartArrayGridOptions = {
    chartId: string;
    cells: ChartArrayCell[];
    getViewId?: (chartId: string, cell: ChartArrayCell, index: number) => string;
};

export type UseChartArrayGridResult = {
    layoutRef: RefObject<ChartArrayLayoutHandle | null>;
    scrollContainerRef: RefObject<HTMLDivElement | null>;
    metrics: ChartArrayMetrics;
    visibleCellIndices: number[];
    viewIds: string[];
    visibleViewIds: string[];
    configuredCellCount: number;
    cellCount: number;
    cellKeys: string[];
    getCellBounds: (index: number) => ChartArrayRect | undefined;
    hasValidVisibleCells: boolean;
    hasCanvas: boolean;
};

export function useChartArrayGrid({
    chartId,
    cells,
    getViewId = (id, cell, index) => getChartArrayViewId(id, cell.key, index),
}: UseChartArrayGridOptions): UseChartArrayGridResult {
    const layoutRef = useRef<ChartArrayLayoutHandle>(null);
    const scrollContainerRef = useRef<HTMLDivElement>(null);
    const cellCount = cells.length;
    const configuredCellCount = cellCount;

    const metrics = useChartArrayMetrics(layoutRef, cellCount);
    const visibleCellIndices = useChartArrayVisibleIndices(scrollContainerRef, layoutRef, cellCount);
    const visibleCellIndicesKey = visibleCellIndices.join("\u0000");
    const stableVisibleCellIndices = useMemo(
        () => visibleCellIndices,
        [visibleCellIndicesKey],
    );

    const viewIds = useMemo(
        () => cells.map((cell, index) => getViewId(chartId, cell, index)),
        [cells, chartId, getViewId],
    );

    const visibleViewIds = useMemo(
        () =>
            stableVisibleCellIndices
                .map((index) => viewIds[index])
                .filter((viewId): viewId is string => Boolean(viewId)),
        [stableVisibleCellIndices, viewIds],
    );

    const cellKeys = useMemo(() => cells.map((cell) => cell.key), [cells]);

    const getCellBounds = (index: number) => metrics.cellBounds[index];

    const hasValidVisibleCells =
        stableVisibleCellIndices.length > 0 &&
        stableVisibleCellIndices.every((index) => {
            const bounds = metrics.cellBounds[index];
            return (
                bounds !== undefined &&
                bounds.width > 0 &&
                bounds.height > 0 &&
                Number.isFinite(bounds.width) &&
                Number.isFinite(bounds.height)
            );
        });

    const hasCanvas =
        metrics.rootSize.width > 0 && metrics.rootSize.height > 0 && hasValidVisibleCells;

    return {
        layoutRef,
        scrollContainerRef,
        metrics,
        visibleCellIndices: stableVisibleCellIndices,
        viewIds,
        visibleViewIds,
        configuredCellCount,
        cellCount,
        cellKeys,
        getCellBounds,
        hasValidVisibleCells,
        hasCanvas,
    };
}
