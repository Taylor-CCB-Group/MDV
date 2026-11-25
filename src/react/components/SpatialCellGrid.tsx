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
import { deckHostLayerId } from "@/react/spatialdata/host_overlay_ids";
import { useChartID, useConfig } from "../hooks";
import { type ChartArrayCell, type UseChartArrayGridResult, useChartArrayGrid } from "../hooks/useChartArrayGrid";
import type { ColorGridCells } from "../hooks/useColorGridCells";
import type { DensityGridCells } from "../hooks/useDensityGridCells";
import type { ViewControllerOptions } from "../hooks/useGateLayers";
import type { LayerFillGridCells } from "../hooks/useLayerFillGridCells";
import type { TooltipContent } from "../hooks/useOuterContainerDeckTooltip";
import { useSpatialLayers } from "../spatial_context";
import DeckCellGrid, { type DeckCellLabel } from "./DeckCellGrid";
import { SHARED_CELL_LAYERS_ABOVE_ID, sharedCellLayersGroup } from "./densityGridSharedLayers";
import { getDensityGridViewId } from "./densityGridUtils";
import { splitRenderStackAt } from "./layerFillGridUtils";
import { type CrosshairConfig, isCrosshairShown } from "./pointerCrosshair";
import { useLayersForDevice } from "./scatterDeckShared";
import { useColorGridCellLayers } from "./useColorGridCellLayers";
import { useDensityGridCellLayers } from "./useDensityGridCellLayers";

type SpatialRenderer = ReturnType<typeof useSpatialCanvasRendererFromLayerInputs>;

/** What the SpatialData chart hands its canvas to draw a grid, besides the cells. */
export type SpatialGridHost = {
    device: DeckDevice;
    deviceProps?: Record<string, unknown>;
    onDeviceInitialized: (device: Device) => void;
    onError: (error: Error, layer?: Layer | null) => void;
    controllerOptions: ViewControllerOptions;
    getTooltipContent: (info: PickingInfo) => TooltipContent | undefined;
    setKeyboardActive: (active: boolean) => void;
    onPointerLeave: () => void;
};

/** The SpatialData chart's grids, by what makes their cells. */
export type SpatialCellGridSpec =
    | { by: "densityFields"; host: SpatialGridHost; cells: DensityGridCells }
    | {
          by: "color_fields";
          host: SpatialGridHost;
          cells: ColorGridCells;
          /** The render stack's points layer, as the adapter resolves it; absent while hidden. */
          pointsLayer: Layer | undefined;
      }
    | {
          by: "layer_fill";
          host: SpatialGridHost;
          /** The cells, and the copies of the layer the renderer draws for them (`withLayerFillCopies`). */
          cells: LayerFillGridCells;
      };

/** The render stack entry of MDV's data points. */
const POINTS_ENTRY_ID = deckHostLayerId("scatter");

/**
 * The view token the image layers carry. It names no cell, as they draw in all of them; with no
 * `viewportId`, their tiles load in whichever cell draws them.
 */
const SHARED_IMAGE_VIEW = "shared";

const getCellViewId = (chartId: string, cell: ChartArrayCell, index: number) =>
    getDensityGridViewId(chartId, cell.key, index);

/**
 * MDV overlays drawn over the shared stack rather than in it: the selection layer listens on each
 * cell's canvas, and the gate labels are dragged and clicked.
 */
const CELL_OVERLAY_ENTRY_IDS = new Set([deckHostLayerId("selection"), deckHostLayerId("gate_labels")]);

/**
 * The stack's layers in the stack's order, with its groups built, as `VivSpatialViewer` places
 * them: the renderer's deck layers and the image layers, each under the entry it belongs to.
 * Layers no entry claims go on top. The `CELL_OVERLAY_ENTRY_IDS` overlays are left out, as
 * `overlays`, in the order the renderer gave them.
 */
export function composeStackLayers(
    deckLayers: readonly Layer[],
    images: Map<string, Layer[]>,
    renderStack: RenderStack | undefined,
    layerOrder: readonly string[],
): { shared: Layer[]; overlays: Layer[] } {
    const { records, overlays } = stackRecords(deckLayers, images);
    return { shared: composeRecords(records, renderStack, layerOrder), overlays };
}

type StackRecord = { layer: Layer; entryId: string };

function stackRecords(deckLayers: readonly Layer[], images: Map<string, Layer[]>) {
    const overlays = deckLayers.filter((layer) => CELL_OVERLAY_ENTRY_IDS.has(layer.id));
    const records: StackRecord[] = [
        ...deckLayers
            .filter((layer) => !CELL_OVERLAY_ENTRY_IDS.has(layer.id))
            .map((layer) => ({ layer, entryId: layer.id })),
        ...[...images].flatMap(([entryId, layers]) => layers.map((layer) => ({ layer, entryId }))),
    ];
    return { records, overlays };
}

/**
 * As `composeStackLayers`, in two parts either side of one entry (see `splitRenderStackAt`), for
 * a grid whose cells each draw their own copy of it in between. Layers no entry claims go on top.
 */
export function composeStackLayersAround(
    deckLayers: readonly Layer[],
    images: Map<string, Layer[]>,
    renderStack: RenderStack | undefined,
    layerOrder: readonly string[],
    entryId: string,
): { below: Layer[]; above: Layer[]; overlays: Layer[] } {
    const { records, overlays } = stackRecords(deckLayers, images);
    const split = renderStack ? splitRenderStackAt(renderStack, entryId) : null;
    if (split) {
        const belowIds = new Set(renderStackOrder(split.below, []));
        return {
            below: composeRecords(
                records.filter((record) => belowIds.has(record.entryId)),
                split.below,
                [],
            ),
            above: composeRecords(
                records.filter((record) => !belowIds.has(record.entryId)),
                split.above,
                [],
            ),
            overlays,
        };
    }
    const at = layerOrder.indexOf(entryId);
    const isBelow = (record: StackRecord) => {
        const index = layerOrder.indexOf(record.entryId);
        return at >= 0 && index >= 0 && index < at;
    };
    return {
        below: composeRecords(records.filter(isBelow), undefined, layerOrder),
        above: composeRecords(
            records.filter((record) => !isBelow(record)),
            undefined,
            layerOrder,
        ),
        overlays,
    };
}

function composeRecords(
    records: StackRecord[],
    renderStack: RenderStack | undefined,
    layerOrder: readonly string[],
): Layer[] {
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

type SpatialCellGridProps = {
    renderer: SpatialRenderer;
    renderStack: RenderStack | undefined;
    viewState: SpatialCanvasViewState;
    onViewStateChange: (next: SpatialCanvasViewState) => void;
    onHover: (info: PickingInfo, canvas: HTMLCanvasElement | undefined) => void;
};

/** The SpatialData chart as a grid of cells, of whichever kind `spec` asks for. */
export default function SpatialCellGrid({ spec, ...props }: SpatialCellGridProps & { spec: SpatialCellGridSpec }) {
    if (spec.by === "color_fields") {
        return <SpatialColorGrid host={spec.host} colorCells={spec.cells} pointsLayer={spec.pointsLayer} {...props} />;
    }
    if (spec.by === "layer_fill") {
        return <SpatialLayerFillGrid host={spec.host} layerCells={spec.cells} {...props} />;
    }
    return <SpatialDensityGrid host={spec.host} densityCells={spec.cells} {...props} />;
}

const NO_CELL_LAYERS: readonly Layer[] = [];

/**
 * A cell per fill colour column (or category) of a shapes or labels layer: the render stack in
 * every cell, but for that layer, of which each cell draws its own copy where the layer was.
 */
const SpatialLayerFillGrid = observer(function SpatialLayerFillGrid({
    host,
    layerCells,
    ...props
}: SpatialCellGridProps & { host: SpatialGridHost; layerCells: LayerFillGridCells }) {
    const chartId = useChartID();
    const grid = useChartArrayGrid({ chartId, cells: layerCells.cells, getViewId: getCellViewId });
    const { cells, copies, layerId } = layerCells;
    const visibleCells = useMemo(
        () => grid.visibleCellIndices.filter((index) => cells[index] && grid.viewIds[index]),
        [grid.visibleCellIndices, grid.viewIds, cells],
    );
    const labels = useMemo(() => cells.map((cell) => ({ text: cell.label, color: cell.color })), [cells]);
    const cellEntry = useMemo(
        () => ({ entryId: layerId, copyIds: new Set(copies.map((copy) => copy.id)) }),
        [layerId, copies],
    );
    return (
        <StackCellGrid
            host={host}
            grid={grid}
            labels={labels}
            visibleCells={visibleCells}
            cellLayers={NO_CELL_LAYERS}
            cellEntry={cellEntry}
            {...props}
        />
    );
});

/**
 * Each density field in its own cell: the render stack in every cell, a cell's density layer only
 * in its own.
 */
const SpatialDensityGrid = observer(function SpatialDensityGrid({
    host,
    densityCells,
    ...props
}: SpatialCellGridProps & { host: SpatialGridHost; densityCells: DensityGridCells }) {
    const chartId = useChartID();
    const grid = useChartArrayGrid({ chartId, cells: densityCells.cells, getViewId: getCellViewId });
    const { visibleCells, densityLayers, labels } = useDensityGridCellLayers(grid, densityCells.densityFields);
    return (
        <StackCellGrid
            host={host}
            grid={grid}
            labels={labels}
            visibleCells={visibleCells}
            cellLayers={densityLayers}
            {...props}
        />
    );
});

/**
 * A cell per colour column, or per category of the first: the render stack in every cell but for
 * its data points, which each cell draws its own way (see `useColorGridCellLayers`).
 */
const SpatialColorGrid = observer(function SpatialColorGrid({
    host,
    colorCells,
    pointsLayer,
    ...props
}: SpatialCellGridProps & { host: SpatialGridHost; colorCells: ColorGridCells; pointsLayer: Layer | undefined }) {
    const chartId = useChartID();
    const { recolorPoints } = useSpatialLayers().scatterProps;
    const grid = useChartArrayGrid({ chartId, cells: colorCells.cells, getViewId: getCellViewId });
    const { visibleCells, pointLayers, labels } = useColorGridCellLayers(grid, colorCells, pointsLayer, recolorPoints);
    return (
        <StackCellGrid
            host={host}
            grid={grid}
            labels={labels}
            visibleCells={visibleCells}
            cellLayers={pointLayers}
            ownPointsEntryId={POINTS_ENTRY_ID}
            {...props}
        />
    );
});

/**
 * The SpatialData chart as small multiples (see `DeckCellGrid`). The images, shapes, points and
 * MDV overlays of the render stack are drawn in every cell, once a frame for all of them; each
 * cell's own layers (ids prefixed with its view id) only in it.
 */
const StackCellGrid = observer(function StackCellGrid({
    host,
    grid,
    labels,
    visibleCells,
    cellLayers: ownCellLayers,
    ownPointsEntryId,
    cellEntry,
    renderer,
    renderStack,
    viewState,
    onViewStateChange,
    onHover,
}: SpatialCellGridProps & {
    host: SpatialGridHost;
    grid: UseChartArrayGridResult;
    labels: readonly DeckCellLabel[];
    visibleCells: readonly number[];
    /** The cells' own layers, each id prefixed with its cell's view id. */
    cellLayers: readonly Layer[];
    /** The stack entry the cells draw their own copies of, left out of what every cell draws. */
    ownPointsEntryId?: string;
    /**
     * A stack entry the renderer draws a copy of per cell (each id under its cell's view id), in
     * its place in the stack: the stack is drawn once below the copies and once above them.
     */
    cellEntry?: { entryId: string; copyIds: ReadonlySet<string> };
}) {
    const config = useConfig<CrosshairConfig>();

    const imageLayers = useMemo(
        () => vivImageLayers(renderer.vivLayerProps, { viewId: SHARED_IMAGE_VIEW }),
        [renderer.vivLayerProps],
    );
    const stackDeckLayers = useMemo(
        () =>
            ownPointsEntryId
                ? renderer.deckLayers.filter((layer) => layer.id !== ownPointsEntryId)
                : renderer.deckLayers,
        [renderer.deckLayers, ownPointsEntryId],
    );
    const { below, copies, above, overlays } = useMemo(() => {
        if (!cellEntry) {
            const composed = composeStackLayers(stackDeckLayers, imageLayers, renderStack, renderer.layerOrder);
            const none: Layer[] = [];
            return { below: composed.shared, copies: none, above: none, overlays: composed.overlays };
        }
        const isCopy = (id: string) => cellEntry.copyIds.has(id);
        const composed = composeStackLayersAround(
            stackDeckLayers.filter((layer) => !isCopy(layer.id)),
            new Map([...imageLayers].filter(([entryId]) => !isCopy(entryId))),
            renderStack,
            renderer.layerOrder,
            cellEntry.entryId,
        );
        const copyLayers = [
            ...stackDeckLayers.filter((layer) => isCopy(layer.id)),
            ...[...imageLayers].flatMap(([entryId, layers]) => (isCopy(entryId) ? layers : [])),
        ];
        return { ...composed, copies: copyLayers };
    }, [stackDeckLayers, imageLayers, renderStack, renderer.layerOrder, cellEntry]);
    // The chart drew these on its single canvas until the cells were ready, and a layer belongs to
    // the first deck that draws it: this deck gets copies.
    const sharedBelow = useLayersForDevice(below, host.device);
    const sharedAbove = useLayersForDevice(above, host.device);
    // each cell's own over the stack (gate outlines included, for now), and the overlays over that
    const cellLayers = useLayersForDevice([...ownCellLayers, ...overlays], host.device);
    const cellCopies = useLayersForDevice(copies, host.device);
    // drawn once a frame for all the cells (see `sharedCellLayersGroup`); a cell's copies of an
    // entry go between the stack below the entry and the stack above it
    const layers = cellEntry
        ? [
              sharedCellLayersGroup(sharedBelow),
              ...cellCopies,
              sharedCellLayersGroup(sharedAbove, SHARED_CELL_LAYERS_ABOVE_ID),
              ...cellLayers,
          ]
        : [sharedCellLayersGroup(sharedBelow), ...cellLayers];

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
