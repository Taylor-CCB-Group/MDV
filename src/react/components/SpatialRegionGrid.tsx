import type { Device } from "@luma.gl/core";
import { viewStateFromBounds } from "@spatialdata/core";
import type { RenderStack } from "@spatialdata/layers";
import {
    type RenderStackLayerInputs,
    resolvePickedEntryId,
    type SpatialFeaturePickEvent,
    useSpatialCanvasRendererFromLayerInputs,
    type VivImagePassthroughOptions,
    vivImageLayers,
} from "@spatialdata/vis";
import type { Layer, OrthographicViewState, PickingInfo } from "deck.gl";
import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useMemo, useRef, useState, useSyncExternalStore } from "react";
import type { DeckDevice } from "@/react/spatialdata/deck_backend";
import { createImageLayerRegistry } from "@/react/spatialdata/image_layer_registry";
import { type AssociatedLayerInputs, useAssociatedShapesLayerInputs } from "@/react/spatialdata/table_association";
import SpatialLayer from "@/webgl/SpatialLayer";
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
    /**
     * A SpatialData feature under the pointer: the event, the canvas its `x` / `y` are relative
     * to, and the MDV tooltip items for the cell's region.
     */
    onFeatureHover: (
        event: SpatialFeaturePickEvent,
        anchor: HTMLElement | undefined,
        tooltipItems: CellTooltipItems,
    ) => void;
    onFeatureHoverEnd: () => void;
};

type CellTooltipItems = AssociatedLayerInputs["mdvTooltipItems"];
/** A cell's hover handler: a pick in its view, with the canvas it was made in. */
type CellHover = (info: PickingInfo, canvas: HTMLCanvasElement | undefined) => boolean;

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
    registerHover,
    onFeatureHover,
    pointsLayer,
    onReadError,
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
    registerHover: (viewId: string, hover: CellHover | null) => void;
    onFeatureHover: SpatialRegionGridHost["onFeatureHover"];
    /** The region's data points, as the render stack's points entry (see `regionPointsLayer`). */
    pointsLayer: Layer | undefined;
    /** Why the region's store could not be read, or null once it can. */
    onReadError: (regionKey: string, message: string | null) => void;
}) {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const reader = chart.spatialDataReader;
    const getSnapshot = useCallback(() => reader.read(cell.url).state, [reader, cell.url]);
    const { spatialData, error } = useSyncExternalStore(reader.subscribe, getSnapshot);
    const errorMessage = error ? error.message || String(error) : null;
    useEffect(() => onReadError(cell.key, errorMessage), [onReadError, cell.key, errorMessage]);
    const { layers, layerOrder, mdvTooltipItems } = useAssociatedShapesLayerInputs(
        spatialData ?? undefined,
        layerInputs,
        cell.file,
    );
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
    // the points go through the stack like the region's own layers, so their group's opacity and
    // blending apply as in the single view
    const cellLayers = useMemo(() => {
        const deckLayers = pointsLayer ? [...renderer.deckLayers, pointsLayer] : renderer.deckLayers;
        const { shared } = composeStackLayers(deckLayers, images, renderStack, renderer.layerOrder);
        return shared.map((layer) => inCell(viewId, layer));
    }, [renderer.deckLayers, pointsLayer, images, renderStack, renderer.layerOrder, viewId]);

    useEffect(() => onLayers(viewId, cellLayers), [onLayers, viewId, cellLayers]);
    useEffect(() => () => onLayers(viewId, null), [onLayers, viewId]);

    // A pick in this cell: its layers' ids carry the cell's view id, which the renderer's don't.
    const hover = useRef<CellHover>(() => false);
    hover.current = (info, canvas) => {
        const pickedId = resolvePickedEntryId(info);
        const entryId = pickedId.startsWith(`${viewId}-`) ? pickedId.slice(viewId.length + 1) : pickedId;
        const event = renderer.getFeaturePickEvent(entryId, { index: info.index, object: info.object });
        if (!event) return false;
        onFeatureHover(
            { ...event, coordinateSystem: cell.coordinateSystem, spatialData, entryId, pickInfo: info },
            canvas,
            mdvTooltipItems,
        );
        return true;
    };
    useEffect(() => {
        registerHover(viewId, (info, canvas) => hover.current(info, canvas));
        return () => registerHover(viewId, null);
    }, [registerHover, viewId]);

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
 * A region's data points: the render stack's points layer (its entry's id, opacity and handlers)
 * with only that region's rows and highlighted rows. Picks name a position in the region's rows;
 * the layer's handlers expect one in all of them, so they are given that.
 */
function regionPointsLayer(points: Layer, { rows, positions }: RegionRows, highlightedData: number[]): Layer {
    const toAllRows = (info: PickingInfo): PickingInfo =>
        info.index < 0 ? info : { ...info, index: positions[info.index] ?? -1 };
    const { onHover, onClick } = points.props;
    return points.clone({
        data: rows,
        // the highlight is drawn apart from `data`, so it has to be narrowed as well
        ...(points instanceof SpatialLayer && { highlightedData }),
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
 * Not yet in cells: selection tools, gates, density and filtered-out (grey) points.
 */
export default observer(function SpatialRegionGrid({
    host,
    layerInputs,
    renderStack,
    vivPassthrough,
    pointsLayer,
}: {
    host: SpatialRegionGridHost;
    /** The render stack's SpatialData layers, as for the single view, without the MDV overlays. */
    layerInputs: RenderStackLayerInputs;
    renderStack: RenderStack | undefined;
    vivPassthrough: VivImagePassthroughOptions;
    /**
     * The render stack's points layer over every region's rows, as the adapter resolves it (with
     * its entry's id and inherited opacity); absent while the points are hidden.
     */
    pointsLayer: Layer | undefined;
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
    // a store that could not be read: the cell says so, and offers to read it again
    const [readErrors, setReadErrors] = useState<Readonly<Record<string, string>>>({});
    const onReadError = useCallback((regionKey: string, message: string | null) => {
        setReadErrors((current) => {
            if (message) return current[regionKey] === message ? current : { ...current, [regionKey]: message };
            if (!(regionKey in current)) return current;
            const { [regionKey]: _removed, ...rest } = current;
            return rest;
        });
    }, []);
    const retry = useCallback(() => chart.spatialDataReader.retry(), [chart]);
    const labels = useMemo<DeckCellLabel[]>(
        () =>
            cells.map((cell) => {
                const count = rowsByRegion?.get(cell.key)?.rows.length;
                const readError = readErrors[cell.key];
                return {
                    text: count === undefined ? cell.label : `${cell.label} (${count})`,
                    color: LABEL_COLOR,
                    ...(readError && {
                        notice: {
                            text: `Could not read this region's SpatialData: ${readError}`,
                            action: { label: "Retry", onClick: retry },
                        },
                    }),
                };
            }),
        [cells, rowsByRegion, readErrors, retry],
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
    // each region's points, with only its own highlighted rows: the highlight is drawn apart from `data`
    const highlighted = pointsLayer instanceof SpatialLayer ? pointsLayer.props.highlightedData : undefined;
    const highlightedByRegion = useMemo(() => {
        if (!Array.isArray(highlighted) || highlighted.length === 0 || !regionColumn?.data || !regionColumn.values) {
            return null;
        }
        return partitionRowsByRegion(highlighted, regionColumn.data, regionColumn.values);
    }, [highlighted, regionColumn]);
    const regionPoints = useMemo(() => {
        const byRegion = new Map<string, Layer>();
        if (!pointsLayer) return byRegion;
        for (const index of visibleCells) {
            const cell = cells[index];
            if (!cell) continue;
            const rows = rowsByRegion?.get(cell.key) ?? { rows: NO_ROWS, positions: NO_ROWS };
            const regionHighlighted = Array.from(highlightedByRegion?.get(cell.key)?.rows ?? []);
            byRegion.set(cell.key, regionPointsLayer(pointsLayer, rows, regionHighlighted));
        }
        return byRegion;
    }, [pointsLayer, visibleCells, cells, rowsByRegion, highlightedByRegion]);

    const layers = useMemo(
        () => visibleCells.flatMap((index) => cellLayers[grid.viewIds[index] ?? ""] ?? []),
        [visibleCells, grid.viewIds, cellLayers],
    );
    const deviceLayers = useLayersForDevice(layers, host.device);

    // each cell handles a hover in its view: only it knows what its layers show
    const hoverHandlers = useRef(new Map<string, CellHover>());
    const registerHover = useCallback((viewId: string, hover: CellHover | null) => {
        if (hover) hoverHandlers.current.set(viewId, hover);
        else hoverHandlers.current.delete(viewId);
    }, []);
    const { onFeatureHoverEnd } = host;
    const onHover = useCallback(
        (info: PickingInfo, canvas: HTMLCanvasElement | undefined) => {
            const hover = info.picked ? hoverHandlers.current.get(info.viewport?.id ?? "") : undefined;
            // deck reports moving off a feature as an unpicked hover
            if (!hover?.(info, canvas)) onFeatureHoverEnd();
        },
        [onFeatureHoverEnd],
    );

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
                        registerHover={registerHover}
                        onFeatureHover={host.onFeatureHover}
                        pointsLayer={regionPoints.get(cell.key)}
                        onReadError={onReadError}
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
                onHover={onHover}
                setKeyboardActive={host.setKeyboardActive}
                onPointerLeave={host.onPointerLeave}
            />
        </>
    );
});
