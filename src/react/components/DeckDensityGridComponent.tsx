import type { Deck, Layer, OrthographicView, OrthographicViewState, PickingInfo } from "@deck.gl/core";
import { action } from "mobx";
import { observer } from "mobx-react-lite";
import { useCallback } from "react";
import { getPickingInfoWithAlternates } from "@/lib/deckPicking";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import { useChart } from "../context";
import type { DualContourLegacyConfig } from "../contour_state";
import { useChartID, useConfig, useFilteredIndices } from "../hooks";
import { type ChartArrayCell, useChartArrayGrid } from "../hooks/useChartArrayGrid";
import type { DensityGridCells } from "../hooks/useDensityGridCells";
import useGateLayers from "../hooks/useGateLayers";
import { useSpatialLayers } from "../spatial_context";
import DeckCellGrid from "./DeckCellGrid";
import type { DeckScatterConfig, DeckScatterReact } from "./DeckScatterReactWrapper";
import { getDensityGridViewId } from "./densityGridUtils";
import { useLayersForDevice, useScatterDeckDevice, useZoomOnFilter } from "./scatterDeckShared";
import { useDensityGridCellLayers } from "./useDensityGridCellLayers";

const getCellViewId = (chartId: string, cell: ChartArrayCell, index: number) =>
    getDensityGridViewId(chartId, cell.key, index);

/**
 * The deck scatter plot with each density field in its own cell (see `DeckCellGrid`). The points,
 * gates and selection are drawn in every cell; a cell's density layer only in its own.
 */
export default observer(function DeckDensityGridComponent({ cells, densityFields }: DensityGridCells) {
    const chartId = useChartID();
    const chart = useChart<DeckScatterConfig, DeckScatterReact>();
    const config = useConfig<DeckScatterConfig & DualContourLegacyConfig>();
    const rows = useFilteredIndices();
    const {
        scatterProps: { scatterplotLayer, greyScatterplotLayer, getTooltip, setScatterKeyboardActive },
        selectionLayer,
    } = useSpatialLayers();
    const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();
    const grid = useChartArrayGrid({ chartId, cells, getViewId: getCellViewId });

    // every cell is the same size, so the first stands for them all
    const cellWidth = grid.metrics.cellBounds[0]?.width ?? 0;
    const cellHeight = grid.metrics.cellBounds[0]?.height ?? 0;
    useZoomOnFilter(rows, cellWidth, cellHeight, { refitOnResize: true });

    const { visibleCells, densityLayers, labels } = useDensityGridCellLayers(grid, densityFields);

    const layers: Layer[] = [
        greyScatterplotLayer,
        scatterplotLayer,
        ...densityLayers,
        gateDisplayLayer,
        selectionLayer,
        gateLabelLayer,
    ].filter((layer): layer is NonNullable<typeof layer> => layer != null);
    const { device, deviceProps, onDeviceInitialized, onError } = useScatterDeckDevice(layers.map((layer) => layer.id));
    // deck picks synchronously only on WebGL; on WebGPU the hover pick is all there is
    const pickAlternates = device.status === "ready" && device.backend === "webgl";
    const deckLayers = useLayersForDevice(layers, device);

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
            dragPan={controllerOptions.dragPan}
            getTooltipContent={getTooltipContent}
            setKeyboardActive={setScatterKeyboardActive}
            exportTarget={chart}
            crosshair={config.grid?.crosshair}
        />
    );
});
