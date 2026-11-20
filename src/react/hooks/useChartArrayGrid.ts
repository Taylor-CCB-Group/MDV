import { type RefObject, useCallback, useMemo, useRef, useState } from "react";
import type { ChartArrayLayoutHandle } from "../components/ChartArrayLayout";
import type { ChartArrayRect } from "../components/chartArrayUtils";
import { getChartArrayViewId } from "../components/chartArrayUtils";
import { type ChartArrayMetrics, useChartArrayMetrics } from "./useChartArrayMetrics";
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
    /** The layout, to read; give `ChartArrayLayout` `setLayout` as its ref. */
    layoutRef: RefObject<ChartArrayLayoutHandle | null>;
    setLayout: (handle: ChartArrayLayoutHandle | null) => void;
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

const defaultGetViewId = (chartId: string, cell: ChartArrayCell, index: number) =>
    getChartArrayViewId(chartId, cell.key, index);

export function useChartArrayGrid({
    chartId,
    cells,
    getViewId = defaultGetViewId,
}: UseChartArrayGridOptions): UseChartArrayGridResult {
    const layoutRef = useRef<ChartArrayLayoutHandle>(null);
    const scrollContainerRef = useRef<HTMLDivElement>(null);
    const cellCount = cells.length;
    const configuredCellCount = cellCount;

    // The layout can be replaced while this hook stays mounted (a grid shows a placeholder
    // while its device changes), so count its mounts: the observers follow the new elements.
    const [layoutMounts, setLayoutMounts] = useState(0);
    const setLayout = useCallback((handle: ChartArrayLayoutHandle | null) => {
        layoutRef.current = handle;
        if (handle) setLayoutMounts((n) => n + 1);
    }, []);

    const cellKeys = useMemo(() => cells.map((cell) => cell.key), [cells]);
    // what the observers watch: a new layout, or a new cell in a slot (a new element even with
    // the same number of cells), means observing again
    const cellsKey = `${layoutMounts}\u0001${cellKeys.join("\u0000")}`;
    const metrics = useChartArrayMetrics(layoutRef, cellCount, cellsKey);
    const visibleCellIndices = useChartArrayVisibleIndices(scrollContainerRef, layoutRef, cellCount, cellsKey);
    const visibleCellIndicesKey = visibleCellIndices.join("\u0000");
    // biome-ignore lint/correctness/useExhaustiveDependencies: a new array only when the indices change
    const stableVisibleCellIndices = useMemo(() => visibleCellIndices, [visibleCellIndicesKey]);

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

    const getCellBounds = useCallback((index: number) => metrics.cellBounds[index], [metrics.cellBounds]);

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

    const hasCanvas = metrics.rootSize.width > 0 && metrics.rootSize.height > 0 && hasValidVisibleCells;

    return {
        layoutRef,
        setLayout,
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
