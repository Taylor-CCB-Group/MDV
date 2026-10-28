import DeckGL, { type DeckGLRef } from "@deck.gl/react";
import { type Layer, OrthographicView, OrbitView, type PickingInfo } from "@deck.gl/core";
import { observer } from "mobx-react-lite";
import { useChartSize, useConfig, useFilteredIndices, useParamColumns } from "../hooks";
import { LineLayer } from "@deck.gl/layers";
import { useCallback, useEffect, useId, useMemo, useRef, useState } from "react";
import { useChart } from "../context";
import type { DeckScatterConfig, DeckScatterReact } from "./DeckScatterReactWrapper";
import { action } from "mobx";
import type { DataColumn, LoadedDataColumn, NumberDataType } from "@/charts/charts";
import { allNumeric } from "@/lib/columnTypeHelpers";
import { SpatialAnnotationProvider, useSpatialLayers } from "../spatial_context";
import SelectionOverlay from "./SelectionOverlay";
import { useScatterRadius } from "../scatter_state";
import AxisComponent from "./AxisComponent";
import { useRebindDeckEvents } from "../hooks/useRebindDeckEvents";
import useGateLayers from "../hooks/useGateLayers";
import FieldContourLegend from "./FieldContourLegend";
import {
    useCategoryContourLegend,
    useFieldContourLegend,
    type DualContourLegacyConfig,
    type HoveredContour,
} from "../contour_state";
import { getPickingInfoWithAlternates } from "@/lib/deckPicking";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import { useOuterContainerDeckTooltip } from "../hooks/useOuterContainerDeckTooltip";
import {
    type DeckDevice,
    type DeckErrorHandlers,
    resolveDeckBackend,
    useDeckDevice,
    useDeckDeviceProps,
} from "@/react/spatialdata/deck_backend";
import { LayerFailureTracker } from "@/react/spatialdata/layer_failures";
import { ErrorBoundary } from "react-error-boundary";
import ErrorComponentReactWrapper from "./ErrorComponentReactWrapper";

const errorMessage = (error: unknown) => (error instanceof Error ? error.message : String(error));

//todo this should be in a common place etc.
const colMid = ({ minMax }: DataColumn<NumberDataType>) => minMax[0] + (minMax[1] - minMax[0]) / 2;

// const margin = { top: 10, right: 10, bottom: 40, left: 60 };
/** todo this should be common for viv / scatter_state, pending refactor
 * there should be hooks getting the range of a filtered column
 * so that multiple charts can re-use the computation (but if not used, it's not computed)
 */
function useZoomOnFilter(data: Uint32Array) {
    const [width, height] = useChartSize();
    const [cx, cy, cz] = useParamColumns();
    // const data = useFilteredIndices(); //<< does calling this multiple times mean wasting space?
    const config = useConfig<DeckScatterConfig>();
    const chart = useChart() as any;
    const { pendingRecenter } = chart;
    useEffect(() => {
        if (chart.ignoreStateUpdate) return;
        if (data.length === 0) return; // [0, 0, 1, 1];
        if (!pendingRecenter && !config.zoom_on_filter) return;
        //there is also cx.minMax, cy.minMax - but not for filtered indices
        let minX = Number.POSITIVE_INFINITY;
        let maxX = Number.NEGATIVE_INFINITY;
        let minY = Number.POSITIVE_INFINITY;
        let maxY = Number.NEGATIVE_INFINITY;
        // todo z... probably refactor for less repetition
        // ^ maybe we could have properties on the columns that are the min/max of the filtered data
        for (let i = 0; i < data.length; i++) {
            const d = data[i];
            const x = cx.data[d];
            const y = cy.data[d];
            if (!Number.isFinite(x) || !Number.isFinite(y)) {
                console.warn("undefined data in scatterplot");
                continue;
            }
            if (x < minX) minX = x;
            if (x > maxX) maxX = x;
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
        }
        action(() => {
            chart.ignoreStateUpdate = true;
            chart.pendingRecenter = false;

            // Calculate range
            const dx = maxX - minX;
            const dy = maxY - minY;

            // Use relative epsilon based on data scale or absolute minimum
            const epsilon = Math.max(1e-9, Math.min(Math.abs(minX), Math.abs(maxX)) * 1e-6);
            const safeDx = Math.max(epsilon, dx);
            const safeDy = Math.max(epsilon, dy);

            // Calculate zoom with fallback and NaN check
            let zoom = Math.log2(Math.min(width / safeDx, height / safeDy)) - 0.6;
            if (!Number.isFinite(zoom)) {
                console.warn("Invalid zoom value detected, using fallback");
                zoom = 0; // Default zoom fallback
            }

            config.viewState = {
                target: [(minX + maxX) / 2, (minY + maxY) / 2, 0],
                zoom,
                // any kind of interpolator would probably mean changing how we handle viewState
                // transitionInterpolator: new FlyToInterpolator({speed: 1.2}),
            };
            if (config.dimension === "3d") {
                const vs = config.viewState;
                vs.rotationOrbit = vs.rotationOrbit || 0;
                vs.rotationX = vs.rotationX || 0;
                if (!config.zoom_on_filter) {
                    if (!cz) {
                        console.warn("3D viewState requires z coordinate");
                        return;
                    }
                    const params = [cx, cy, cz];
                    if (!allNumeric(params)) {
                        throw new Error("3D viewState requires numeric coordinates");
                    }
                    vs.target = params.map(colMid) as [number, number, number];
                }
            }

            chart.ignoreStateUpdate = false;
        })();
    }, [data, cx, cy, cz, width, height, config, pendingRecenter, chart]);
}

/**
 * The device deck draws with (see `deck_backend.ts`), and what to do when it fails. On
 * WebGPU, a chart following the default backend moves to WebGL when a layer, the device or
 * the canvas fails; one set to WebGPU explicitly stays on it, with the errors in the console
 * (each kind once: a WebGPU validation error otherwise repeats every frame).
 */
function useScatterDeckDevice(layerIds: readonly string[]) {
    const chart = useChart<DeckScatterConfig, DeckScatterReact>();
    const config = useConfig<DeckScatterConfig>();
    const { device, fallBack } = useDeckDevice(resolveDeckBackend(config.render_backend));
    const activeBackend = device.status === "ready" ? device.backend : null;
    const fallbackReason = device.status === "ready" ? device.fallbackReason : undefined;
    useEffect(() => {
        chart.setActiveBackend(activeBackend, fallbackReason);
        return () => chart.setActiveBackend(null);
    }, [chart, activeBackend, fallbackReason]);
    const layerIdsRef = useRef(layerIds);
    layerIdsRef.current = layerIds;
    /** Fall back to WebGL, if this chart is on WebGPU because it follows the default. */
    const failedOnWebGPU = useCallback(
        (reason: string) => {
            // after the switch, the old WebGPU device's errors are stale
            if (chart.activeBackend !== "webgpu" || chart.config.render_backend != null) return;
            console.error(`[MDV] ${reason}; using WebGL`);
            fallBack(reason);
        },
        [chart, fallBack],
    );
    const failureTracker = useMemo(
        () =>
            new LayerFailureTracker({
                entryIds: () => layerIdsRef.current,
                fail: (layerId, failure) => failedOnWebGPU(`${layerId} failed on WebGPU: ${failure.reason}`),
                canvas: (message) => failedOnWebGPU(`the canvas failed on WebGPU: ${message}`),
            }),
        [failedOnWebGPU],
    );
    useEffect(() => () => failureTracker.dispose(), [failureTracker]);
    const handlers = useMemo<DeckErrorHandlers>(
        () => ({
            device: (error) => failureTracker.device(error),
            layer: (error, layerId) => failureTracker.layer(error, layerId),
        }),
        [failureTracker],
    );
    const deckProps = useDeckDeviceProps(device, fallBack, handlers);
    return { device, failedOnWebGPU, ...deckProps };
}

/**
 * A deck layer instance belongs to the first Deck that draws it: its models and buffers are
 * on that Deck's device. These layers are memoised outside the canvas, so a Deck remounted
 * for another backend would inherit the old device's state and fail to draw. Give each
 * device its own copies, one per layer instance.
 */
function useLayersForDevice(layers: Layer[], device: DeckDevice): Layer[] {
    // biome-ignore lint/correctness/useExhaustiveDependencies: a new device is the trigger
    const copies = useMemo(() => new WeakMap<Layer, Layer>(), [device]);
    return layers.map((layer) => {
        const cached = copies.get(layer);
        if (cached) return cached;
        const copy: Layer = layer.clone({});
        copies.set(layer, copy);
        return copy;
    });
}

/**
 * Currently this implementation is somewhat separate from VivMDVReact / scatter_state,
 * but it should be more unified in the future.
 *
 * Things that need work:
 * - GeoJsonEditableLayer & associated selection state
 *  - think about selection in 3d
 * - density (in future, more different types of layers)
 *  - fancy ways of e.g. selecting things based on density (metaballs)
 * - fix mouse events in popout
 * - auto-scaling of radius (including when parameter changes)
 * - different scaling along different axes when not using physical units
 * - highlighting
 * - tooltips
 *   - enhanced version with more info (in general for all charts, or at least new ones)
 * - axis configuration (including via direct manipulation)
 */
const DeckScatter = observer(function DeckScatterComponent({
    setHoveredContour,
}: {
    setHoveredContour: (contour: HoveredContour | null) => void;
}) {
    const id = useId();
    const [width, height] = useChartSize();
    const [cx, cy, ...density] = useParamColumns() as LoadedDataColumn<"double">[];
    const cz = useParamColumns()[2] as LoadedDataColumn<"double">;
    const data = useFilteredIndices(); //changed to fallback to simplerFilteredIndices when filterColumn is not set
    const config = useConfig<DeckScatterConfig>();
    const contourConfig = useConfig<DualContourLegacyConfig>();
    const { viewState, dimension } = config;
    const is2d = dimension === "2d";
    //todo more clarity on radius units - but large radius was causing big problems after deck upgrade
    const radiusScale = useScatterRadius();
    //todo colorBy should be done differently (also bearing in mind multiple layers)

    //todo this shouldn't be repeated here and in AxisComponent
    const xSize = is2d ? config.axis.x.size : 0;
    const ySize = is2d ? config.axis.y.size : 0;

    const margin = useMemo(
        () =>
            is2d
                ? {
                      top: 10,
                      right: 10,
                      bottom: xSize,
                      left: ySize,
                  }
                : {
                      top: 0,
                      right: 0,
                      bottom: 0,
                      left: 0,
                  },
        [is2d, xSize, ySize],
    );
    // never negative (a chart narrower than its axes): WebGPU rejects a negative viewport
    const chartWidth = Math.max(0, width - margin.left - margin.right);
    //there could be a potential off-by-one/two error somewhere down the line
    //if we don't fully understand reasons for `- 3.5` here.
    //prevents overlapping with x-axis.
    const chartHeight = Math.max(0, height - margin.top - margin.bottom - 3.5);
    useZoomOnFilter(data);

    const { scatterProps, selectionLayer } = useSpatialLayers();
    // this is now somewhat able to render for "2d", pending further tweaks
    //! beware unproject from here is not what we want, should review
    const { scatterplotLayer, greyScatterplotLayer, getTooltip, setScatterKeyboardActive } = scatterProps;

    const {
        gateLabelLayer,
        gateDisplayLayer,
        controllerOptions,
    } = useGateLayers();

    const legendFields = useFieldContourLegend(contourConfig.densityFields);
    const legendCategories = useCategoryContourLegend(contourConfig);
    const showLegend = contourConfig.field_legend.display;
    const legendPosition = { x: 10, y: 10 };

    const axisLinesLayer = useMemo(() => {
        if (is2d) return null;
        return new LineLayer({
            id: `axis-lines-${id}`,
            data: [
                { sourcePosition: [0, 0, 0], targetPosition: [cx.minMax[1], 0, 0], color: [255, 0, 0] },
                { sourcePosition: [0, 0, 0], targetPosition: [0, cy.minMax[1], 0], color: [0, 255, 0] },
                { sourcePosition: [0, 0, 0], targetPosition: [0, 0, cz.minMax[1]], color: [0, 0, 255] },
            ],
            getSourcePosition: (d: any) => d.sourcePosition,
            getTargetPosition: (d: any) => d.targetPosition,
            getColor: (d: any) => d.color,
            getWidth: 1,
            updateTriggers: {
                getSourcePosition: [cx.minMax, cy.minMax, cz.minMax],
                getTargetPosition: [cx.minMax, cy.minMax, cz.minMax],
            },
        });
    }, [is2d, id, cx.minMax, cy.minMax, cz?.minMax]);

    // we need an OrthographicView to prevent wrapping etc...
    // if in future we have subgraphs sharing a canvas, we will need to
    // make sure that the view is set up correctly for each subgraph.
    const view = useMemo(() => {
        return config.dimension === "2d" ? new OrthographicView({
            id: `scatterplot-view-${id}`,
            // controller: true,
            width: chartWidth,
            height: chartHeight,
            x: 0,
            y: 0,
            flipY: false,
        }) : new OrbitView({
            id: `scatterplot-view-${id}`,
            // controller: true,
            width: chartWidth,
            height: chartHeight,
            x: 0,
            y: 0,
        });
    }, [chartWidth, chartHeight, config.dimension, id]);

    //! deck doesn't like it if we change the layers array - better to toggle visibility
    const layers = [greyScatterplotLayer, scatterplotLayer, gateDisplayLayer, selectionLayer, gateLabelLayer, axisLinesLayer,
    ].filter((layer): layer is NonNullable<typeof layer> => layer != null);
    const { device, failedOnWebGPU, deviceProps, onDeviceInitialized, onError } = useScatterDeckDevice(
        layers.map((layer) => layer.id),
    );
    // deck picks synchronously only on WebGL; on WebGPU the hover pick is all there is
    const pickAlternates = device.status === "ready" && device.backend === "webgl";
    const deckLayers = useLayersForDevice(layers, device);

    const deckContainerRef = useRef<HTMLDivElement | null>(null);
    const deckRef = useRef<DeckGLRef | null>(null);
    const getTooltipContent = useCallback(
        (info: PickingInfo) => {
            const richInfo = getPickingInfoWithAlternates(info, pickAlternates ? deckRef.current?.deck : null);
            return getCombinedScatterTooltip(
                richInfo,
                {
                    gateDisplayLayerId: gateDisplayLayer?.id,
                    gateLabelLayerId: gateLabelLayer?.id,
                    getPointTooltip: getTooltip,
                },
            );
        },
        [gateDisplayLayer?.id, gateLabelLayer?.id, getTooltip, pickAlternates],
    );
    const {
        clearTooltip,
        getTooltip: getPortalTooltip,
        suppressTooltipUntilPointerUp,
        tooltipPortal,
    } = useOuterContainerDeckTooltip(getTooltipContent, deckContainerRef);

    // unproject used for updating ranges - use deck viewport instead of layer
    const unproject = useCallback((coords: [number, number]) => {
        if (!deckRef.current?.deck) {
            throw new Error("Deck instance not yet initialized");
        }
        const deck = deckRef.current.deck;
        // we may want to deal with multiple viewports for "splatter-plot" & other scenarios.
        const viewport = deck.getViewports()[0];
        if (!viewport) {
            throw new Error("No viewport available");
        }
        // Unproject from screen coordinates to world coordinates
        return viewport.unproject([coords[0], coords[1]]);
    }, []);
    useRebindDeckEvents(deckRef);

    // we want default controller options, but we want a new one when the outerContainer changes
    // this doesn't seem to help re-register mouse events.
    // const controller = useMemo(() => ({inertia: 10+Math.random()}), [outerContainer])


    return (
        <>
            <AxisComponent config={config} unproject={unproject}>
                <div
                    ref={deckContainerRef}
                    aria-label="Scatter plot"
                    style={{ width: "100%", height: "100%", outline: "none" }}
                    onPointerDown={suppressTooltipUntilPointerUp}
                    onMouseDown={() => {
                        setScatterKeyboardActive(true);
                    }}
                    onMouseEnter={() => setScatterKeyboardActive(true)}
                    onMouseLeave={() => {
                        clearTooltip();
                        setScatterKeyboardActive(false);
                    }}
                >
                    {device.status === "pending" ? (
                        <div className="h-full w-full p-2">Starting WebGPU…</div>
                    ) : (
                        <ErrorBoundary
                            // deck picks its device once, when created: remount on a backend change.
                            key={device.backend}
                            onError={(error) => {
                                if (device.backend === "webgpu") failedOnWebGPU(`the canvas failed on WebGPU: ${errorMessage(error)}`);
                            }}
                            fallbackRender={({ error }) => (
                                <ErrorComponentReactWrapper
                                    error={{
                                        message: errorMessage(error),
                                        ...(error instanceof Error && error.stack ? { stack: error.stack } : {}),
                                    }}
                                    title="Error showing scatter plot canvas"
                                />
                            )}
                        >
                            <DeckGL
                                ref={deckRef}
                                layers={deckLayers}
                                {...(deviceProps && { deviceProps })}
                                onDeviceInitialized={onDeviceInitialized}
                                onError={onError}
                                useDevicePixels={true}
                                controller={{
                                    // double-click finishes polygons / removes vertices in selection tools
                                    doubleClickZoom: false,
                                    dragPan: controllerOptions.dragPan,
                                }}
                                viewState={viewState}
                                // initialViewState={viewState} //consider not using react state for this        
                                views={view}
                                onViewStateChange={(v) => {
                                    action(() => (config.viewState = v.viewState))();
                                }}
                                getTooltip={getPortalTooltip}
                                getCursor={({ isDragging }) => {
                                    return isDragging ? "grabbing" : "crosshair";
                                }}
                            />
                        </ErrorBoundary>
                    )}
                </div>
            </AxisComponent>
            {tooltipPortal}
            {showLegend && (
                <FieldContourLegend
                    fields={legendFields}
                    categories={legendCategories.items}
                    categoryLabel={legendCategories.label}
                    position={legendPosition}
                    onContourHover={setHoveredContour}
                />
            )}
        </>
    );
});

export default () => {
    const chart = useChart();
    const [hoveredContour, setHoveredContour] = useState<HoveredContour | null>(null);
    // in order for SelectionOverlay to work, we need to review how our implementation works
    // vs useScatterplotLayer in spatial_context.tsx
    return (
        <SpatialAnnotationProvider chart={chart} hoveredContour={hoveredContour} setHoveredContour={setHoveredContour}>
            <SelectionOverlay />
            <DeckScatter setHoveredContour={setHoveredContour} />
        </SpatialAnnotationProvider>
    );
};
