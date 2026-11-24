import type { Device } from "@luma.gl/core";
import { viewStateFromBounds } from "@spatialdata/core";
import type { RenderStack } from "@spatialdata/layers";
import {
    type RenderStackLayerInputs,
    useSpatialCanvasRendererFromLayerInputs,
    type VivImagePassthroughOptions,
    vivImageLayers,
} from "@spatialdata/vis";
import type { Layer, OrthographicViewState, PickingInfo } from "deck.gl";
import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useMemo, useState, useSyncExternalStore } from "react";
import type { DeckDevice } from "@/react/spatialdata/deck_backend";
import { createImageLayerRegistry } from "@/react/spatialdata/image_layer_registry";
import { useAssociatedShapesLayerInputs } from "@/react/spatialdata/table_association";
import { useChart } from "../context";
import { useChartID, useConfig, useFieldSpec } from "../hooks";
import { type ChartArrayCell, useChartArrayGrid } from "../hooks/useChartArrayGrid";
import type { ViewControllerOptions } from "../hooks/useGateLayers";
import type { TooltipContent } from "../hooks/useOuterContainerDeckTooltip";
import DeckCellGrid, { type DeckCellLabel } from "./DeckCellGrid";
import { getDensityGridViewId } from "./densityGridUtils";
import {
    moveRegionViews,
    partitionRowsByRegion,
    type RegionRows,
    type RegionViewState,
    type SpatialRegionCell,
    spatialRegionCells,
    toRegionViewState,
} from "./regionGridUtils";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "./SpatialDataMDVReact";
import { composeStackLayers } from "./SpatialDensityGrid";
import { useLayersForDevice } from "./scatterDeckShared";

/** What the SpatialData chart hands its region grid, besides the layers. */
export type SpatialRegionGridHost = {
    device: DeckDevice;
    deviceProps?: Record<string, unknown>;
    onDeviceInitialized: (device: Device) => void;
    onError: (error: Error, layer?: Layer | null) => void;
    controllerOptions: ViewControllerOptions;
    getTooltipContent: (info: PickingInfo) => TooltipContent | undefined;
    setKeyboardActive: (active: boolean) => void;
    onPointerLeave: () => void;
};

const LABEL_COLOR = [226, 232, 240] as const;
const NO_ROWS = new Uint32Array(0);

// the deck cell grid draws a layer only in the cell whose view id starts its id
const getCellViewId = (chartId: string, cell: ChartArrayCell, index: number) =>
    getDensityGridViewId(chartId, cell.key, index);

/**
 * A region's view as deck's, the same object for the same state: a cell given an unchanged view
 * does not render again when another cell moves.
 */
const orthographicViews = new WeakMap<RegionViewState, OrthographicViewState>();
function toOrthographic(state: RegionViewState): OrthographicViewState {
    let view = orthographicViews.get(state);
    if (!view) {
        view = { target: [state.target[0], state.target[1], 0], zoom: state.zoom };
        orthographicViews.set(state, view);
    }
    return view;
}

/** Gives a layer an id under a cell's view id, so it draws only in that cell. */
function inCell(viewId: string, layer: Layer): Layer {
    return layer.clone({ id: `${viewId}-${layer.id}` });
}

/**
 * The layers of one region's cell: its SpatialData drawn with the chart's layers, each under the
 * cell's view id. Draws nothing itself; it reports the layers (and a view that fits them, until
 * the cell has one) to the grid, which draws every cell from one deck. Memoised (as an observer):
 * the grid renders again whenever a cell reports, and a cell rendered again with the same props
 * would report again.
 */
const RegionCellLayers = observer(function RegionCellLayers({
    cell,
    viewId,
    layerInputs,
    renderStack,
    vivPassthrough,
    width,
    height,
    viewState,
    onLayers,
    onFit,
    publishesImageData,
}: {
    cell: SpatialRegionCell;
    viewId: string;
    layerInputs: RenderStackLayerInputs;
    renderStack: RenderStack | undefined;
    vivPassthrough: VivImagePassthroughOptions;
    width: number;
    height: number;
    viewState: OrthographicViewState | undefined;
    onLayers: (viewId: string, layers: Layer[] | null) => void;
    onFit: (regionKey: string, viewState: OrthographicViewState) => void;
    /** Whether the layer dialog's image controls (histograms) read this cell's images. */
    publishesImageData: boolean;
}) {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const reader = chart.spatialDataReader;
    const getSnapshot = useCallback(() => reader.read(cell.url).state, [reader, cell.url]);
    const { spatialData } = useSyncExternalStore(reader.subscribe, getSnapshot);
    const { layers, layerOrder } = useAssociatedShapesLayerInputs(spatialData ?? undefined, layerInputs, cell.file);
    const spatialViewState = useMemo(
        () =>
            viewState
                ? {
                      target: [viewState.target?.[0] ?? 0, viewState.target?.[1] ?? 0] as [number, number],
                      zoom: typeof viewState.zoom === "number" ? viewState.zoom : 0,
                  }
                : null,
        [viewState],
    );
    const renderer = useSpatialCanvasRendererFromLayerInputs({
        spatialData,
        coordinateSystem: cell.coordinateSystem,
        layerInputs: { layers, layerOrder },
        renderOrder: layerOrder,
        viewState: spatialViewState,
        width,
        height,
        sortDeckLayers: true,
        autoFit: false,
        vivPassthrough,
    });
    // tiles follow this cell's view, which no other cell shares
    const images = useMemo(
        () => vivImageLayers(renderer.vivLayerProps, { viewId, viewportId: viewId }),
        [renderer.vivLayerProps, viewId],
    );
    const cellLayers = useMemo(() => {
        const { shared } = composeStackLayers(renderer.deckLayers, images, renderStack, renderer.layerOrder);
        return shared.map((layer) => inCell(viewId, layer));
    }, [renderer.deckLayers, images, renderStack, renderer.layerOrder, viewId]);

    useEffect(() => onLayers(viewId, cellLayers), [onLayers, viewId, cellLayers]);
    useEffect(() => () => onLayers(viewId, null), [onLayers, viewId]);

    // The layer dialog's image controls read one store's images, as in the single view. The stack
    // is edited in place, so a new generation means a new registry, as there.
    const stack = chart.config.renderStack;
    const generation = chart.renderStackGeneration;
    const { getImageLoadedDataByElementKey, getLayerLoadState } = renderer;
    useEffect(() => {
        void generation;
        if (!publishesImageData || !stack) return;
        const registry = createImageLayerRegistry(stack, getImageLoadedDataByElementKey, getLayerLoadState, cell.key);
        chart.setImageLayerRegistry(registry);
        return () => {
            // another cell may have published since
            if (chart.imageLayerRegistry === registry) chart.setImageLayerRegistry(undefined);
        };
    }, [chart, publishesImageData, stack, generation, getImageLoadedDataByElementKey, getLayerLoadState, cell.key]);

    // fit the cell to its region once its layers have bounds, as the single view does
    const { getWorldBoundsForVisibleLayers, hasEnabledLayers, hasLayersDrawn, isBlocking } = renderer;
    useEffect(() => {
        if (viewState || !hasEnabledLayers || isBlocking || !hasLayersDrawn) return;
        if (width <= 0 || height <= 0) return;
        const bounds = getWorldBoundsForVisibleLayers();
        if (!bounds) return;
        const fitted = viewStateFromBounds(bounds, width, height);
        onFit(cell.key, { target: [fitted.target[0], fitted.target[1], 0], zoom: fitted.zoom });
    }, [
        viewState,
        hasEnabledLayers,
        isBlocking,
        hasLayersDrawn,
        width,
        height,
        getWorldBoundsForVisibleLayers,
        onFit,
        cell.key,
    ]);
    return null;
});

/**
 * A region's data points: the chart's points layer with only that region's rows, under the cell's
 * view id. Picks name a position in the region's rows; the layer's handlers expect one in all of
 * them, so they are given that.
 */
function regionPointsLayer(points: Layer, viewId: string, { rows, positions }: RegionRows): Layer {
    const toAllRows = (info: PickingInfo): PickingInfo =>
        info.index < 0 ? info : { ...info, index: positions[info.index] ?? -1 };
    const { onHover, onClick } = points.props;
    return points.clone({
        id: `${viewId}-${points.id}`,
        data: rows,
        onHover: onHover
            ? (info: PickingInfo, event: Parameters<typeof onHover>[1]) => onHover(toAllRows(info), event)
            : null,
        onClick: onClick
            ? (info: PickingInfo, event: Parameters<typeof onClick>[1]) => onClick(toAllRows(info), event)
            : null,
    });
}

/**
 * The SpatialData chart with a cell per region. Every cell draws its region's SpatialData with the
 * chart's render stack, as a template (each layer names an element, found in each region's store),
 * and the region's data points. Each cell is framed on its own region and pans and zooms alone.
 * Only the cells in view load and draw anything.
 *
 * Not yet in cells: selection tools, gates, density, filtered-out (grey) points and feature tooltips.
 */
export default observer(function SpatialRegionGrid({
    host,
    layerInputs,
    renderStack,
    vivPassthrough,
    pointsLayer,
    pointsVisible,
}: {
    host: SpatialRegionGridHost;
    /** The render stack's SpatialData layers, as for the single view, without the MDV overlays. */
    layerInputs: RenderStackLayerInputs;
    renderStack: RenderStack | undefined;
    vivPassthrough: VivImagePassthroughOptions;
    /** The chart's points layer over every region's rows. */
    pointsLayer: Layer | undefined;
    pointsVisible: boolean;
}) {
    const chartId = useChartID();
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const config = useConfig<SpatialDataMdvReactConfig>();
    const regions = chart.dataStore.regions;
    const regionColumn = useFieldSpec(regions?.region_field);

    const pointRows = pointsLayer?.props.data;
    const rowsByRegion = useMemo(() => {
        if (!(pointRows instanceof Uint32Array) || !regionColumn?.data || !regionColumn.values) return null;
        return partitionRowsByRegion(pointRows, regionColumn.data, regionColumn.values);
    }, [pointRows, regionColumn]);

    const allCells = useMemo(() => spatialRegionCells(regions?.all_regions), [regions]);
    const hideEmpty = config.grid?.hide_empty === true;
    const cells = useMemo(
        () => (hideEmpty && rowsByRegion ? allCells.filter((cell) => rowsByRegion.has(cell.key)) : allCells),
        [allCells, hideEmpty, rowsByRegion],
    );
    const grid = useChartArrayGrid({ chartId, cells, getViewId: getCellViewId });
    const labels = useMemo<DeckCellLabel[]>(
        () =>
            cells.map((cell) => {
                const count = rowsByRegion?.get(cell.key)?.rows.length;
                return { text: count === undefined ? cell.label : `${cell.label} (${count})`, color: LABEL_COLOR };
            }),
        [cells, rowsByRegion],
    );

    // Each cell's view, by region, on the chart (see `regionViewStates`), so it outlives the grid
    // and is saved with the view.
    const { regionViewStates } = chart;
    const linkViews = config.grid?.link_views === true;
    const onFit = useCallback(
        (regionKey: string, viewState: OrthographicViewState) => {
            if (chart.regionViewStates[regionKey]) return;
            chart.setRegionViewStates({ ...chart.regionViewStates, [regionKey]: toRegionViewState(viewState) });
        },
        [chart],
    );
    const viewStates = useMemo(
        () =>
            Object.fromEntries(
                cells.flatMap((cell, index) => {
                    const state = regionViewStates[cell.key];
                    const viewId = grid.viewIds[index];
                    return state && viewId ? [[viewId, toOrthographic(state)]] : [];
                }),
            ),
        [cells, regionViewStates, grid.viewIds],
    );
    const onViewStateChange = useCallback(
        (next: OrthographicViewState, viewId: string) => {
            const cell = cells[grid.viewIds.indexOf(viewId)];
            if (!cell) return;
            chart.setRegionViewStates(moveRegionViews(regionViewStates, cell.key, toRegionViewState(next), linkViews));
        },
        [chart, cells, grid.viewIds, regionViewStates, linkViews],
    );

    const [cellLayers, setCellLayers] = useState<Readonly<Record<string, Layer[]>>>({});
    const onLayers = useCallback((viewId: string, layers: Layer[] | null) => {
        setCellLayers((current) => {
            if (layers) return current[viewId] === layers ? current : { ...current, [viewId]: layers };
            if (!(viewId in current)) return current;
            const { [viewId]: _removed, ...rest } = current;
            return rest;
        });
    }, []);

    const visibleCells = grid.visibleCellIndices;
    // the image controls read the chart's own region when it is in view, else the first cell in view
    const imageDataCell = visibleCells.find((index) => cells[index]?.key === config.region) ?? visibleCells[0];
    const layers = useMemo(
        () =>
            visibleCells.flatMap((index) => {
                const viewId = grid.viewIds[index];
                const cell = cells[index];
                if (!viewId || !cell) return [];
                const regionRows = rowsByRegion?.get(cell.key);
                const points =
                    pointsVisible && pointsLayer
                        ? [regionPointsLayer(pointsLayer, viewId, regionRows ?? { rows: NO_ROWS, positions: NO_ROWS })]
                        : [];
                return [...(cellLayers[viewId] ?? []), ...points];
            }),
        [visibleCells, grid.viewIds, cells, rowsByRegion, pointsVisible, pointsLayer, cellLayers],
    );
    const deviceLayers = useLayersForDevice(layers, host.device);

    const fallbackViewState = useMemo<OrthographicViewState>(() => ({ target: [0, 0, 0], zoom: 0 }), []);
    const { getTooltipContent } = host;
    const tooltipContent = useCallback((info: PickingInfo) => getTooltipContent(info), [getTooltipContent]);

    return (
        <>
            {visibleCells.map((index) => {
                const cell = cells[index];
                const viewId = grid.viewIds[index];
                const bounds = grid.metrics.cellBounds[index];
                if (!cell || !viewId || !bounds) return null;
                return (
                    <RegionCellLayers
                        key={viewId}
                        cell={cell}
                        viewId={viewId}
                        layerInputs={layerInputs}
                        renderStack={renderStack}
                        vivPassthrough={vivPassthrough}
                        width={bounds.width}
                        height={bounds.height}
                        viewState={viewStates[viewId]}
                        onLayers={onLayers}
                        onFit={onFit}
                        publishesImageData={index === imageDataCell}
                    />
                );
            })}
            <DeckCellGrid
                grid={grid}
                labels={labels}
                visibleCells={visibleCells}
                layers={deviceLayers}
                device={host.device}
                deviceProps={host.deviceProps}
                onDeviceInitialized={host.onDeviceInitialized}
                onError={host.onError}
                viewState={fallbackViewState}
                viewStates={viewStates}
                onViewStateChange={onViewStateChange}
                controllerOptions={host.controllerOptions}
                getTooltipContent={tooltipContent}
                setKeyboardActive={host.setKeyboardActive}
                onPointerLeave={host.onPointerLeave}
            />
        </>
    );
});
