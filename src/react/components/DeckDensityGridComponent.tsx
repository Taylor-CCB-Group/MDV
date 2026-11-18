import { Deck, type Layer, OrthographicView, type OrthographicViewState, type PickingInfo } from "@deck.gl/core";
import type { Device } from "@luma.gl/core";
import { action } from "mobx";
import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useMemo, useRef } from "react";
import { syncEditableLayerEvents } from "@/lib/deckMonkeypatch";
import { getPickingInfoWithAlternates } from "@/lib/deckPicking";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import { DensityOverlayLayer } from "@/webgl/SpatialLayer";
import { useChart } from "../context";
import { type DualContourLegacyConfig, useFieldContour } from "../contour_state";
import { getFieldColor } from "../fieldColorManager";
import { useChartID, useConfig, useFilteredIndices } from "../hooks";
import { type ChartArrayCell, useChartArrayGrid } from "../hooks/useChartArrayGrid";
import { useDensityGridCells } from "../hooks/useDensityGridCells";
import useGateLayers from "../hooks/useGateLayers";
import { useOuterContainerDeckTooltip } from "../hooks/useOuterContainerDeckTooltip";
import { useRebindDeckEvents } from "../hooks/useRebindDeckEvents";
import { useRange, useSpatialLayers } from "../spatial_context";
import ChartArrayLayout from "./ChartArrayLayout";
import type { DeckScatterConfig, DeckScatterReact } from "./DeckScatterReactWrapper";
import { getDensityGridViewId, shouldDrawInDensityGridCell } from "./densityGridUtils";
import { useLayersForDevice, useScatterDeckDevice, useZoomOnFilter } from "./scatterDeckShared";

const getCellViewId = (chartId: string, cell: ChartArrayCell, index: number) =>
    getDensityGridViewId(chartId, cell.key, index);

const centred = (message: string) => <div className="flex h-full items-center justify-center text-sm">{message}</div>;

type ExportCell = { element: HTMLElement; canvas: HTMLCanvasElement | undefined; text: string; color: number[] };

/**
 * The grid as one image: each cell's canvas where the cell is, with its field name, as the
 * labels on screen are HTML. Cells scrolled out of view have no drawing, only their label.
 */
function composeGridImage(root: HTMLElement, cells: readonly ExportCell[]) {
    const scale = globalThis.devicePixelRatio || 1;
    const rootRect = root.getBoundingClientRect();
    const out = document.createElement("canvas");
    out.width = Math.round(rootRect.width * scale);
    out.height = Math.round(rootRect.height * scale);
    const ctx = out.getContext("2d");
    if (!ctx) return null;
    ctx.scale(scale, scale);
    ctx.font = "12px sans-serif";
    ctx.textBaseline = "middle";
    for (const { element, canvas, text, color } of cells) {
        const r = element.getBoundingClientRect();
        const x = r.left - rootRect.left;
        const y = r.top - rootRect.top;
        if (canvas && canvas.width > 0 && canvas.height > 0) ctx.drawImage(canvas, x, y, r.width, r.height);
        const width = Math.min(ctx.measureText(text).width + 16, r.width - 16);
        ctx.fillStyle = "rgba(15, 23, 42, 0.76)";
        ctx.fillRect(x + 8, y + 8, width, 24);
        ctx.fillStyle = `rgb(${color[0]}, ${color[1]}, ${color[2]})`;
        ctx.fillText(text, x + 16, y + 20, width - 16);
    }
    return out;
}

/**
 * The deck scatter plot with each density field in its own cell. Every cell has its own
 * canvas, so the browser lays them out and scrolls them; one deck, on one device, draws into
 * the canvases of the cells in view (deck's experimental multi-canvas mode, `_canvases`), with
 * a view per canvas under one shared view state, so panning or zooming any cell moves them
 * all. The points, gates and selection are drawn in every cell; a cell's density layer only
 * in its own (see `shouldDrawInDensityGridCell`).
 */
export default observer(function DeckDensityGridComponent() {
    const chartId = useChartID();
    const chart = useChart<DeckScatterConfig, DeckScatterReact>();
    const config = useConfig<DeckScatterConfig & DualContourLegacyConfig>();
    const rows = useFilteredIndices();
    const { modelMatrix } = useRange();
    const {
        scatterProps: { scatterplotLayer, greyScatterplotLayer, getTooltip, setScatterKeyboardActive },
        selectionLayer,
    } = useSpatialLayers();
    const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();
    const { cells, densityFields, configuredFieldCount } = useDensityGridCells();
    const grid = useChartArrayGrid({ chartId, cells, getViewId: getCellViewId });

    // every cell is the same size, so the first stands for them all
    const cellWidth = grid.metrics.cellBounds[0]?.width ?? 0;
    const cellHeight = grid.metrics.cellBounds[0]?.height ?? 0;
    useZoomOnFilter(rows, cellWidth, cellHeight, { refitOnResize: true });

    const visibleCells = useMemo(
        () => grid.visibleCellIndices.filter((index) => densityFields[index] && grid.viewIds[index]),
        [grid.visibleCellIndices, grid.viewIds, densityFields],
    );
    const visibleFields = useMemo(
        () => visibleCells.map((index) => densityFields[index]),
        [visibleCells, densityFields],
    );
    // aligned with `visibleFields`
    const contours = useFieldContour({
        id: "densityGrid",
        fields: visibleFields,
        fill: config.contour_fill,
        bandwidth: config.contour_bandwidth ?? 0.1,
        intensity: config.contour_intensity ?? 0.1,
        opacity: config.contour_opacity ?? 0.2,
        fillThreshold: config.contour_fillThreshold ?? 2,
    });
    const densityLayers = useMemo(
        () =>
            visibleCells.map(
                (index, i) =>
                    new DensityOverlayLayer({
                        id: `${grid.viewIds[index]}-density`,
                        contourLayers: [contours[i]],
                        modelMatrix,
                    }),
            ),
        [visibleCells, grid.viewIds, contours, modelMatrix],
    );

    // a view per cell in view, drawn into that cell's canvas (the canvas id is the view id)
    const views = useMemo(
        () =>
            visibleCells.map(
                (index) =>
                    new OrthographicView({
                        id: grid.viewIds[index],
                        canvasId: grid.viewIds[index],
                        flipY: false,
                        // double-click finishes polygons / removes vertices in selection tools
                        controller: { doubleClickZoom: false, dragPan: controllerOptions.dragPan },
                    }),
            ),
        [visibleCells, grid.viewIds, controllerOptions.dragPan],
    );
    const viewState = useMemo(
        () => Object.fromEntries(views.map((view) => [view.id, config.viewState])),
        [views, config.viewState],
    );

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

    /** Each cell's canvas, by view id. */
    const canvases = useRef(new Map<string, HTMLCanvasElement>());
    const deckHolder = useRef<{ deck: Deck<OrthographicView[]> | null }>({ deck: null });
    useRebindDeckEvents(deckHolder);

    // deck reads these once, when it is created; route them to the current handlers
    const deviceCallbacks = useRef({ onDeviceInitialized, onError });
    deviceCallbacks.current = { onDeviceInitialized, onError };
    const backend = device.status === "ready" ? device.backend : null;
    // biome-ignore lint/correctness/useExhaustiveDependencies: a new device is the trigger; deviceProps belong to it
    useEffect(() => {
        if (device.status !== "ready") return;
        const deck = new Deck<OrthographicView[]>({
            ...(deviceProps && { deviceProps }),
            onDeviceInitialized: (created: Device) => deviceCallbacks.current.onDeviceInitialized(created),
            onError: (error, layer) => deviceCallbacks.current.onError(error, layer),
            // the cells' canvases are supplied as they come into view
            _canvases: [],
            parent: grid.scrollContainerRef.current,
            useDevicePixels: true,
            // the editable selection layer listens on deck's canvases, which exist from here on
            onLoad: () => syncEditableLayerEvents(deck),
        });
        deckHolder.current.deck = deck;
        return () => {
            deckHolder.current.deck = null;
            deck.finalize();
        };
    }, [device]);

    const getTooltipContent = useCallback(
        (info: PickingInfo) => {
            const richInfo = getPickingInfoWithAlternates(info, pickAlternates ? deckHolder.current.deck : null);
            return getCombinedScatterTooltip(richInfo, {
                gateDisplayLayerId: gateDisplayLayer?.id,
                gateLabelLayerId: gateLabelLayer?.id,
                getPointTooltip: getTooltip,
            });
        },
        [gateDisplayLayer?.id, gateLabelLayer?.id, getTooltip, pickAlternates],
    );
    // a pick's x / y are relative to the canvas it was made in
    const tooltipAnchor = useCallback((info: PickingInfo) => canvases.current.get(info.viewport?.id ?? ""), []);
    const {
        clearTooltip,
        getTooltip: getPortalTooltip,
        suppressTooltipUntilPointerUp,
        tooltipPortal,
    } = useOuterContainerDeckTooltip(getTooltipContent, tooltipAnchor);

    // after every render: the canvases in view, their views, and the layers
    useEffect(() => {
        const deck = deckHolder.current.deck;
        if (!deck) return;
        const drawn = views.filter((view) => canvases.current.has(view.id));
        deck.setProps({
            _canvases: drawn.flatMap((view) => canvases.current.get(view.id) ?? []),
            views: drawn,
            viewState,
            layers: deckLayers,
            layerFilter: ({ layer, viewport }) => shouldDrawInDensityGridCell(layer.id, viewport.id),
            onViewStateChange: ({ viewState: next }: { viewState: OrthographicViewState }) => {
                action(() => {
                    config.viewState = next;
                })();
            },
            getTooltip: getPortalTooltip,
            getCursor: ({ isDragging }) => (isDragging ? "grabbing" : "crosshair"),
        });
        // deck makes an event manager per canvas; follow them as cells scroll in and out
        syncEditableLayerEvents(deck);
    });

    useEffect(() => {
        chart.getPlotCanvas = () => {
            const deck = deckHolder.current.deck;
            const root = grid.layoutRef.current?.root;
            if (!deck || !root) return null;
            // drawing now lets the canvases be read before the browser clears them
            deck.redraw("image export");
            // WebGPU work is submitted at the end of a frame, which this is not part of
            deck.getCanvasContext()?.device?.submit();
            const cellElements = grid.layoutRef.current?.cells ?? [];
            return composeGridImage(
                root,
                densityFields.flatMap((field, index) => {
                    const element = cellElements[index];
                    if (!element) return [];
                    const canvas = canvases.current.get(grid.viewIds[index]);
                    return [{ element, canvas, text: field.name || field.field, color: getFieldColor(field.field) }];
                }),
            );
        };
        return () => {
            chart.getPlotCanvas = undefined;
        };
    }, [chart, grid.layoutRef, grid.viewIds, densityFields]);

    const renderCell = useCallback(
        (index: number) => {
            const field = densityFields[index];
            const viewId = grid.viewIds[index];
            if (!field || !viewId) return null;
            const [r, g, b] = getFieldColor(field.field);
            return (
                <>
                    <canvas
                        // a canvas holds one kind of context for good: a new one for another backend
                        key={backend ?? "pending"}
                        id={viewId}
                        ref={(canvas) => {
                            if (canvas) canvases.current.set(viewId, canvas);
                            else canvases.current.delete(viewId);
                        }}
                    />
                    <div
                        className="pointer-events-none absolute top-2 left-2 z-10 max-w-[calc(100%-16px)] truncate rounded-sm px-2 py-1 text-xs"
                        style={{ backgroundColor: "rgba(15, 23, 42, 0.76)", color: `rgb(${r}, ${g}, ${b})` }}
                        title={field.name || field.field}
                    >
                        {field.name || field.field}
                    </div>
                </>
            );
        },
        [densityFields, grid.viewIds, backend],
    );

    if (configuredFieldCount === 0) return centred("Choose density fields to build the grid.");
    if (densityFields.length === 0) return centred("Loading density fields…");
    if (device.status === "pending") return centred("Starting WebGPU…");

    return (
        <>
            {/* biome-ignore lint/a11y/noStaticElementInteractions: these handlers only track keyboard focus and the tooltip for the canvases inside, as the single plot's do */}
            <div
                ref={grid.scrollContainerRef}
                className="relative h-full w-full min-w-0 overflow-auto overflow-x-hidden"
                onPointerDown={suppressTooltipUntilPointerUp}
                onMouseDown={() => setScatterKeyboardActive(true)}
                onMouseEnter={() => setScatterKeyboardActive(true)}
                onMouseLeave={() => {
                    clearTooltip();
                    setScatterKeyboardActive(false);
                }}
            >
                <ChartArrayLayout
                    cellCount={grid.cellCount}
                    cellKeys={grid.cellKeys}
                    layoutRef={grid.layoutRef}
                    renderCell={renderCell}
                />
            </div>
            {tooltipPortal}
        </>
    );
});
