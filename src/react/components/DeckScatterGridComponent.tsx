import type { Deck, Layer, OrthographicView, OrthographicViewState, PickingInfo } from "@deck.gl/core";
import { action } from "mobx";
import { observer } from "mobx-react-lite";
import { useCallback, useEffect } from "react";
import { getPickingInfoWithAlternates } from "@/lib/deckPicking";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import { useChart } from "../context";
import { useChartID, useConfig, useFilteredIndices } from "../hooks";
import { type ChartArrayCell, type UseChartArrayGridResult, useChartArrayGrid } from "../hooks/useChartArrayGrid";
import type { ColorGridCells } from "../hooks/useColorGridCells";
import type { DensityGridCells } from "../hooks/useDensityGridCells";
import useGateLayers from "../hooks/useGateLayers";
import { useSpatialLayers } from "../spatial_context";
import DeckCellGrid, { type DeckCellLabel } from "./DeckCellGrid";
import type { DeckScatterConfig, DeckScatterReact } from "./DeckScatterReactWrapper";
import { sharedCellLayersGroup } from "./densityGridSharedLayers";
import { getDensityGridViewId } from "./densityGridUtils";
import { isCrosshairShown } from "./pointerCrosshair";
import { useLayersForDevice, useScatterDeckDevice, useZoomOnFilter } from "./scatterDeckShared";
import { useColorGridCellLayers } from "./useColorGridCellLayers";
import { useDensityGridCellLayers } from "./useDensityGridCellLayers";

const withoutNull = (layers: (Layer | null | undefined)[]) => layers.filter((layer): layer is Layer => layer != null);

const getCellViewId = (chartId: string, cell: ChartArrayCell, index: number) =>
    getDensityGridViewId(chartId, cell.key, index);

/**
 * The deck scatter plot with each density field in its own cell (see `DeckCellGrid`). The points,
 * gates and selection are drawn in every cell; a cell's density layer only in its own.
 */
export const DeckDensityGridComponent = observer(function DeckDensityGridComponent({
    cells,
    densityFields,
}: DensityGridCells) {
    const chartId = useChartID();
    const grid = useChartArrayGrid({ chartId, cells, getViewId: getCellViewId });
    const { visibleCells, densityLayers, labels } = useDensityGridCellLayers(grid, densityFields);
    return (
        <ScatterCellGrid
            grid={grid}
            labels={labels}
            visibleCells={visibleCells}
            cellLayers={densityLayers}
            pointsInEveryCell
        />
    );
});

/**
 * The deck scatter plot with a cell per colour column, or per category of the first (see
 * `DeckCellGrid`). Each cell draws its own points; gates, selection and filtered-out (grey) points
 * are drawn in every cell.
 */
export const DeckColorGridComponent = observer(function DeckColorGridComponent({
    colorCells,
}: {
    colorCells: ColorGridCells;
}) {
    const chartId = useChartID();
    const { scatterplotLayer, recolorPoints } = useSpatialLayers().scatterProps;
    const grid = useChartArrayGrid({ chartId, cells: colorCells.cells, getViewId: getCellViewId });
    const { visibleCells, pointLayers, labels } = useColorGridCellLayers(
        grid,
        colorCells,
        scatterplotLayer,
        recolorPoints,
    );
    return (
        <ScatterCellGrid
            grid={grid}
            labels={labels}
            visibleCells={visibleCells}
            cellLayers={pointLayers}
            pointsInEveryCell={false}
        />
    );
});

/**
 * Small multiples of the deck scatter plot. The gates, filtered-out (grey) points and, unless the
 * cells draw their own, the points are drawn once a frame for all the cells (see
 * `sharedCellLayersGroup`); each cell's own layers (ids prefixed with its view id) only in it.
 */
const ScatterCellGrid = observer(function ScatterCellGrid({
    grid,
    labels,
    visibleCells,
    cellLayers: ownCellLayers,
    pointsInEveryCell,
}: {
    grid: UseChartArrayGridResult;
    labels: readonly DeckCellLabel[];
    visibleCells: readonly number[];
    /** The cells' own layers, each id prefixed with its cell's view id. */
    cellLayers: readonly Layer[];
    /** Draw the chart's points in every cell; false when each cell draws its own. */
    pointsInEveryCell: boolean;
}) {
    const chart = useChart<DeckScatterConfig, DeckScatterReact>();
    const config = useConfig<DeckScatterConfig>();
    const rows = useFilteredIndices();
    const {
        scatterProps: { scatterplotLayer, greyScatterplotLayer, getTooltip, setScatterKeyboardActive },
        selectionLayer,
    } = useSpatialLayers();
    const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();

    // turned on just now (not restored with a view): fit the cells, which are smaller than the
    // single plot whose view the chart still has; once measured, `useZoomOnFilter` does it
    useEffect(() => {
        if (!chart.fitGridOnShow) return;
        chart.fitGridOnShow = false;
        action(() => {
            chart.pendingRecenter = true;
        })();
    }, [chart]);
    // every cell is the same size, so the first stands for them all
    const cellWidth = grid.metrics.cellBounds[0]?.width ?? 0;
    const cellHeight = grid.metrics.cellBounds[0]?.height ?? 0;
    useZoomOnFilter(rows, cellWidth, cellHeight, { refitOnResize: true });

    const sharedLayers = withoutNull([
        greyScatterplotLayer,
        pointsInEveryCell ? scatterplotLayer : null,
        gateDisplayLayer,
    ]);
    // over the shared layers: each cell's own, the selection (it listens on each cell's canvas), and
    // the gate labels (dragged and clicked)
    const cellLayers = withoutNull([...ownCellLayers, selectionLayer, gateLabelLayer]);
    const { device, deviceProps, onDeviceInitialized, onError } = useScatterDeckDevice(
        [...sharedLayers, ...cellLayers].map((layer) => layer.id),
    );
    // deck picks synchronously only on WebGL; on WebGPU the hover pick is all there is
    const pickAlternates = device.status === "ready" && device.backend === "webgl";
    const sharedDeckLayers = useLayersForDevice(sharedLayers, device);
    const cellDeckLayers = useLayersForDevice(cellLayers, device);
    const deckLayers = [sharedCellLayersGroup(sharedDeckLayers), ...cellDeckLayers];

    const getTooltipContent = useCallback(
        (info: PickingInfo, deck: Deck<OrthographicView[]> | null) => {
            const richInfo = getPickingInfoWithAlternates(info, pickAlternates ? deck : null);
            return getCombinedScatterTooltip(richInfo, {
                gateDisplayLayerId: gateDisplayLayer?.id,
                gateLabelLayerId: gateLabelLayer?.id,
                getPointTooltip: getTooltip,
            });
        },
        [gateDisplayLayer?.id, gateLabelLayer?.id, getTooltip, pickAlternates],
    );
    const onViewStateChange = useCallback(
        (next: OrthographicViewState) => {
            action(() => {
                config.viewState = next;
            })();
        },
        [config],
    );

    return (
        <DeckCellGrid
            grid={grid}
            labels={labels}
            visibleCells={visibleCells}
            layers={deckLayers}
            device={device}
            deviceProps={deviceProps}
            onDeviceInitialized={onDeviceInitialized}
            onError={onError}
            viewState={config.viewState}
            onViewStateChange={onViewStateChange}
            controllerOptions={controllerOptions}
            getTooltipContent={getTooltipContent}
            setKeyboardActive={setScatterKeyboardActive}
            exportTarget={chart}
            crosshair={isCrosshairShown(config)}
        />
    );
});
