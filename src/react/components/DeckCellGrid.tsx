import { Deck, type Layer, OrthographicView, type OrthographicViewState, type PickingInfo } from "@deck.gl/core";
import type { Device } from "@luma.gl/core";
import { applyWebGPUPickingFix, applyWebGPUPresentationAlphaFix } from "@spatialdata/layers";
import { type ReactNode, useCallback, useEffect, useLayoutEffect, useMemo, useRef } from "react";
import { syncEditableLayerEvents } from "@/lib/deckMonkeypatch";
import type { DeckDevice } from "@/react/spatialdata/deck_backend";
import type { UseChartArrayGridResult } from "../hooks/useChartArrayGrid";
import { type TooltipContent, useOuterContainerDeckTooltip } from "../hooks/useOuterContainerDeckTooltip";
import { useRebindDeckEvents } from "../hooks/useRebindDeckEvents";
import ChartArrayLayout from "./ChartArrayLayout";
import { CrosshairLines } from "./CrosshairArea";
import { shouldDrawInDensityGridCell } from "./densityGridUtils";
import { hideCrosshair, placeCrosshair } from "./pointerCrosshair";

export type DeckCellLabel = { text: string; color: readonly number[] };

type ExportCell = { element: HTMLElement; canvas: HTMLCanvasElement | undefined } & DeckCellLabel;

/**
 * The grid as one image: each cell's canvas where the cell is, with its label, as the labels on
 * screen are HTML. Cells scrolled out of view have no drawing, only their label.
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

export type DeckCellGridProps = {
    grid: UseChartArrayGridResult;
    /** One per cell, by index. */
    labels: readonly DeckCellLabel[];
    /** The cells to draw: those in view. */
    visibleCells: readonly number[];
    /**
     * Every layer, for every cell. A layer whose id starts with a cell's view id draws only in
     * that cell; the rest draw in all of them (see `shouldDrawInDensityGridCell`).
     */
    layers: Layer[];
    device: DeckDevice;
    deviceProps?: Record<string, unknown>;
    onDeviceInitialized: (device: Device) => void;
    onError: (error: Error, layer?: Layer | null) => void;
    /** Shared by every cell. */
    viewState: OrthographicViewState;
    onViewStateChange: (next: OrthographicViewState) => void;
    dragPan: boolean;
    /** The tooltip for a pick; `deck` is there to pick again (for overlapping layers). */
    getTooltipContent: (info: PickingInfo, deck: Deck<OrthographicView[]> | null) => TooltipContent | undefined;
    /** A hover, with the canvas it was made in (`info.x` / `info.y` are relative to it). */
    onHover?: (info: PickingInfo, canvas: HTMLCanvasElement | undefined) => void;
    setKeyboardActive: (active: boolean) => void;
    onPointerLeave?: () => void;
    /** Given `getPlotCanvas`, for exporting the grid as an image. */
    exportTarget?: { getPlotCanvas?: () => HTMLCanvasElement | null };
    /** Mark the pointer's position in every cell with a horizontal and a vertical line. */
    crosshair?: boolean;
};

/**
 * A cell's canvas, made with this window's `document` even when the chart is popped out. deck and
 * luma recognise a canvas with `instanceof HTMLCanvasElement`, against this window's class, and a
 * canvas React makes in a popout's document is the popout window's: luma takes it for something
 * else and never sizes it, so it keeps the default 300 × 150 drawing buffer and its view takes
 * input only in that corner.
 */
function CellCanvas({
    id,
    onCanvas,
}: {
    id: string;
    onCanvas: (id: string, canvas: HTMLCanvasElement | null) => void;
}) {
    const hostRef = useRef<HTMLDivElement>(null);
    useLayoutEffect(() => {
        const host = hostRef.current;
        if (!host) return;
        const canvas = document.createElement("canvas");
        canvas.id = id;
        host.appendChild(canvas);
        onCanvas(id, canvas);
        return () => {
            onCanvas(id, null);
            canvas.remove();
        };
    }, [id, onCanvas]);
    return <div ref={hostRef} className="mdv-chart-array__canvas" />;
}

/**
 * Small multiples from one deck. Every cell has its own canvas, so the browser lays them out and
 * scrolls them; one deck, on one device, draws into the canvases of the cells in view (deck's
 * experimental multi-canvas mode, `_canvases`), with a view per canvas under one shared view
 * state, so panning or zooming any cell moves them all.
 */
export default function DeckCellGrid({
    grid,
    labels,
    visibleCells,
    layers,
    device,
    deviceProps,
    onDeviceInitialized,
    onError,
    viewState,
    onViewStateChange,
    dragPan,
    getTooltipContent,
    onHover,
    setKeyboardActive,
    onPointerLeave,
    exportTarget,
    crosshair = false,
}: DeckCellGridProps) {
    // a view per cell in view, drawn into that cell's canvas (the canvas id is the view id)
    const views = useMemo(
        () =>
            visibleCells.flatMap((index) => {
                const id = grid.viewIds[index];
                if (!id) return [];
                // double-click finishes polygons / removes vertices in selection tools
                return [
                    new OrthographicView({
                        id,
                        canvasId: id,
                        flipY: false,
                        controller: { doubleClickZoom: false, dragPan },
                    }),
                ];
            }),
        [visibleCells, grid.viewIds, dragPan],
    );
    const viewStates = useMemo(() => Object.fromEntries(views.map((view) => [view.id, viewState])), [views, viewState]);

    /** Each cell's canvas, by view id. */
    const canvases = useRef(new Map<string, HTMLCanvasElement>());
    const onCanvas = useCallback((id: string, canvas: HTMLCanvasElement | null) => {
        if (canvas) canvases.current.set(id, canvas);
        else canvases.current.delete(id);
    }, []);
    const deckHolder = useRef<{ deck: Deck<OrthographicView[]> | null }>({ deck: null });
    useRebindDeckEvents(deckHolder);

    // deck reads these once, when it is created; route them to the current handlers
    const callbacks = useRef({ onDeviceInitialized, onError, onViewStateChange, onHover });
    callbacks.current = { onDeviceInitialized, onError, onViewStateChange, onHover };
    const backend = device.status === "ready" ? device.backend : null;
    // biome-ignore lint/correctness/useExhaustiveDependencies: a new device is the trigger; deviceProps belong to it
    useEffect(() => {
        if (device.status !== "ready") return;
        const deck = new Deck<OrthographicView[]>({
            ...(deviceProps && { deviceProps }),
            onDeviceInitialized: (created: Device) => {
                // a canvas per cell: picks there are mirrored, and the canvases opaque, without these
                applyWebGPUPickingFix(created);
                applyWebGPUPresentationAlphaFix(created);
                callbacks.current.onDeviceInitialized(created);
            },
            onError: (error, layer) => callbacks.current.onError(error, layer),
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

    const tooltipContent = useCallback(
        (info: PickingInfo) => getTooltipContent(info, deckHolder.current.deck),
        [getTooltipContent],
    );
    // a pick's x / y are relative to the canvas it was made in
    const tooltipAnchor = useCallback((info: PickingInfo) => canvases.current.get(info.viewport?.id ?? ""), []);
    const {
        clearTooltip,
        getTooltip: getPortalTooltip,
        suppressTooltipUntilPointerUp,
        tooltipPortal,
    } = useOuterContainerDeckTooltip(tooltipContent, tooltipAnchor);

    // after every render: the canvases in view, their views, and the layers
    useEffect(() => {
        const deck = deckHolder.current.deck;
        if (!deck) return;
        const drawn = views.filter((view) => canvases.current.has(view.id));
        deck.setProps({
            _canvases: drawn.flatMap((view) => canvases.current.get(view.id) ?? []),
            views: drawn,
            viewState: viewStates,
            layers,
            layerFilter: ({ layer, viewport }) => shouldDrawInDensityGridCell(layer.id, viewport.id),
            onViewStateChange: ({ viewState: next }: { viewState: OrthographicViewState }) =>
                callbacks.current.onViewStateChange(next),
            onHover: (info) => callbacks.current.onHover?.(info, canvases.current.get(info.viewport?.id ?? "")),
            getTooltip: getPortalTooltip,
            getCursor: ({ isDragging }) => (isDragging ? "grabbing" : "crosshair"),
        });
        // deck makes an event manager per canvas; follow them as cells scroll in and out
        syncEditableLayerEvents(deck);
    });

    useEffect(() => {
        if (!exportTarget) return;
        exportTarget.getPlotCanvas = () => {
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
                labels.flatMap((label, index) => {
                    const element = cellElements[index];
                    if (!element) return [];
                    return [{ element, canvas: canvases.current.get(grid.viewIds[index]), ...label }];
                }),
            );
        };
        return () => {
            exportTarget.getPlotCanvas = undefined;
        };
    }, [exportTarget, grid.layoutRef, grid.viewIds, labels]);

    const renderCell = useCallback(
        (index: number): ReactNode => {
            const label = labels[index];
            const viewId = grid.viewIds[index];
            if (!label || !viewId) return null;
            const [r, g, b] = label.color;
            return (
                <>
                    <CellCanvas
                        // a canvas holds one kind of context for good: a new one for another backend
                        key={backend ?? "pending"}
                        id={viewId}
                        onCanvas={onCanvas}
                    />
                    <div
                        className="pointer-events-none absolute top-2 left-2 z-10 max-w-[calc(100%-16px)] truncate rounded-sm px-2 py-1 text-xs"
                        style={{ backgroundColor: "rgba(15, 23, 42, 0.76)", color: `rgb(${r}, ${g}, ${b})` }}
                        title={label.text}
                    >
                        {label.text}
                    </div>
                    {crosshair && <CrosshairLines />}
                </>
            );
        },
        [labels, grid.viewIds, backend, crosshair, onCanvas],
    );

    if (device.status === "pending") {
        return <div className="flex h-full items-center justify-center text-sm">Starting WebGPU…</div>;
    }

    return (
        <>
            {/* biome-ignore lint/a11y/noStaticElementInteractions: these handlers only track keyboard focus and the tooltip for the canvases inside, as the single plot's do */}
            <div
                ref={grid.scrollContainerRef}
                className="relative h-full w-full min-w-0 overflow-auto overflow-x-hidden"
                onPointerDown={suppressTooltipUntilPointerUp}
                onPointerMove={
                    crosshair
                        ? (event) =>
                              placeCrosshair(
                                  event.currentTarget,
                                  [...canvases.current.values()].find((canvas) => canvas === event.target),
                                  event.clientX,
                                  event.clientY,
                              )
                        : undefined
                }
                onMouseDown={() => setKeyboardActive(true)}
                onMouseEnter={() => setKeyboardActive(true)}
                onMouseLeave={(event) => {
                    hideCrosshair(event.currentTarget);
                    clearTooltip();
                    setKeyboardActive(false);
                    onPointerLeave?.();
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
}
