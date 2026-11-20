import type { Device } from "@luma.gl/core";
import type { RenderStack } from "@spatialdata/layers";
import {
    renderStackOrder,
    renderStackToDeckLayers,
    type ViewState as SpatialCanvasViewState,
    type useSpatialCanvasRendererFromLayerInputs,
    vivImageLayers,
} from "@spatialdata/vis";
import type { Layer, OrthographicViewState, PickingInfo } from "deck.gl";
import { observer } from "mobx-react-lite";
import { useCallback, useMemo } from "react";
import type { DeckDevice } from "@/react/spatialdata/deck_backend";
import type { DualContourLegacyConfig } from "../contour_state";
import { useChartID, useConfig } from "../hooks";
import { type ChartArrayCell, useChartArrayGrid } from "../hooks/useChartArrayGrid";
import type { DensityGridCells } from "../hooks/useDensityGridCells";
import type { ViewControllerOptions } from "../hooks/useGateLayers";
import type { TooltipContent } from "../hooks/useOuterContainerDeckTooltip";
import DeckCellGrid from "./DeckCellGrid";
import { getDensityGridViewId } from "./densityGridUtils";
import { type CrosshairConfig, isCrosshairShown } from "./pointerCrosshair";
import { useLayersForDevice } from "./scatterDeckShared";
import { useDensityGridCellLayers } from "./useDensityGridCellLayers";

type SpatialRenderer = ReturnType<typeof useSpatialCanvasRendererFromLayerInputs>;

/** What the SpatialData chart hands its canvas to draw a density grid. */
export type SpatialDensityGridHost = DensityGridCells & {
    device: DeckDevice;
    deviceProps?: Record<string, unknown>;
    onDeviceInitialized: (device: Device) => void;
    onError: (error: Error, layer?: Layer | null) => void;
    controllerOptions: ViewControllerOptions;
    getTooltipContent: (info: PickingInfo) => TooltipContent | undefined;
    setKeyboardActive: (active: boolean) => void;
    onPointerLeave: () => void;
};

/**
 * The view token the image layers carry. It names no cell, as they draw in all of them; with no
 * `viewportId`, their tiles load in whichever cell draws them.
 */
const SHARED_IMAGE_VIEW = "shared";

const getCellViewId = (chartId: string, cell: ChartArrayCell, index: number) =>
    getDensityGridViewId(chartId, cell.key, index);

/**
 * The stack's layers in the stack's order, with its groups built, as `VivSpatialViewer` places
 * them: the renderer's deck layers and the image layers, each under the entry it belongs to.
 * Layers no entry claims go on top.
 */
function composeStackLayers(
    deckLayers: readonly Layer[],
    images: Map<string, Layer[]>,
    renderStack: RenderStack | undefined,
    layerOrder: readonly string[],
): Layer[] {
    const records = [
        ...deckLayers.map((layer) => ({ layer, entryId: layer.id })),
        ...[...images].flatMap(([entryId, layers]) => layers.map((layer) => ({ layer, entryId }))),
    ];
    if (renderStack) {
        const stackIds = new Set(renderStackOrder(renderStack, []));
        const entryLayers = new Map<string, Layer[]>();
        const unplaced: Layer[] = [];
        for (const { layer, entryId } of records) {
            if (!stackIds.has(entryId)) {
                unplaced.push(layer);
                continue;
            }
            entryLayers.set(entryId, [...(entryLayers.get(entryId) ?? []), layer]);
        }
        return [...renderStackToDeckLayers(renderStack, { entryLayers }), ...unplaced];
    }
    const order = new Map(layerOrder.map((id, index) => [id, index]));
    return records
        .map((record, index) => ({ ...record, index, order: order.get(record.entryId) ?? Number.POSITIVE_INFINITY }))
        .sort((a, b) => a.order - b.order || a.index - b.index)
        .map((record) => record.layer);
}

/**
 * The SpatialData chart with each density field in its own cell (see `DeckCellGrid`). The images,
 * shapes, points and MDV overlays of the render stack are drawn in every cell; a cell's density
 * layer only in its own.
 */
export default observer(function SpatialDensityGrid({
    host,
    renderer,
    renderStack,
    viewState,
    onViewStateChange,
    onHover,
}: {
    host: SpatialDensityGridHost;
    renderer: SpatialRenderer;
    renderStack: RenderStack | undefined;
    viewState: SpatialCanvasViewState;
    onViewStateChange: (next: SpatialCanvasViewState) => void;
    onHover: (info: PickingInfo, canvas: HTMLCanvasElement | undefined) => void;
}) {
    const chartId = useChartID();
    const config = useConfig<DualContourLegacyConfig & CrosshairConfig>();
    const grid = useChartArrayGrid({ chartId, cells: host.cells, getViewId: getCellViewId });
    const { visibleCells, densityLayers, labels } = useDensityGridCellLayers(grid, host.densityFields);

    const imageLayers = useMemo(
        () => vivImageLayers(renderer.vivLayerProps, { viewId: SHARED_IMAGE_VIEW }),
        [renderer.vivLayerProps],
    );
    const stackLayers = useMemo(
        () => composeStackLayers(renderer.deckLayers, imageLayers, renderStack, renderer.layerOrder),
        [renderer.deckLayers, imageLayers, renderStack, renderer.layerOrder],
    );
    // the density cells go over the stack, gates and selection included, for now
    const allLayers = useMemo(() => [...stackLayers, ...densityLayers], [stackLayers, densityLayers]);
    // The chart drew these on its single canvas until the fields loaded, and a layer belongs to
    // the first deck that draws it: this deck gets copies.
    const layers = useLayersForDevice(allLayers, host.device);

    const deckViewState = useMemo<OrthographicViewState>(
        () => ({ target: [viewState.target[0], viewState.target[1], 0], zoom: viewState.zoom }),
        [viewState],
    );
    const onDeckViewStateChange = useCallback(
        (next: OrthographicViewState) => {
            const [x, y] = next.target ?? viewState.target;
            const zoom = typeof next.zoom === "number" ? next.zoom : viewState.zoom;
            onViewStateChange({ target: [x, y], zoom });
        },
        [onViewStateChange, viewState],
    );
    const { getTooltipContent } = host;
    const tooltipContent = useCallback((info: PickingInfo) => getTooltipContent(info), [getTooltipContent]);

    return (
        <DeckCellGrid
            grid={grid}
            labels={labels}
            visibleCells={visibleCells}
            layers={layers}
            device={host.device}
            deviceProps={host.deviceProps}
            onDeviceInitialized={host.onDeviceInitialized}
            onError={host.onError}
            viewState={deckViewState}
            onViewStateChange={onDeckViewStateChange}
            controllerOptions={host.controllerOptions}
            getTooltipContent={tooltipContent}
            onHover={onHover}
            setKeyboardActive={host.setKeyboardActive}
            onPointerLeave={host.onPointerLeave}
            crosshair={isCrosshairShown(config)}
        />
    );
});
