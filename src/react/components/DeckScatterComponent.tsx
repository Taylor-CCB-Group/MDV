import { OrbitView, OrthographicView, type PickingInfo } from "@deck.gl/core";
import { LineLayer } from "@deck.gl/layers";
import DeckGL, { type DeckGLRef } from "@deck.gl/react";
import { action } from "mobx";
import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useId, useMemo, useRef, useState } from "react";
import type { LoadedDataColumn } from "@/charts/charts";
import { getPickingInfoWithAlternates } from "@/lib/deckPicking";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import { useChart } from "../context";
import {
    type DualContourLegacyConfig,
    type HoveredContour,
    useCategoryContourLegend,
    useFieldContourLegend,
} from "../contour_state";
import { useChartSize, useConfig, useFilteredIndices, useParamColumns } from "../hooks";
import useGateLayers from "../hooks/useGateLayers";
import { useOuterContainerDeckTooltip } from "../hooks/useOuterContainerDeckTooltip";
import { useRebindDeckEvents } from "../hooks/useRebindDeckEvents";
import { useScatterRadius } from "../scatter_state";
import { SpatialAnnotationProvider, useSpatialLayers } from "../spatial_context";
import AxisComponent from "./AxisComponent";
import DeckDensityGridComponent from "./DeckDensityGridComponent";
import type { DeckScatterConfig, DeckScatterReact } from "./DeckScatterReactWrapper";
import { isDensityGridShown } from "./densityGridUtils";
import FieldContourLegend from "./FieldContourLegend";
import SelectionOverlay from "./SelectionOverlay";
import { ScatterDeckBoundary, useLayersForDevice, useScatterDeckDevice, useZoomOnFilter } from "./scatterDeckShared";

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
    // the area the points are drawn in, which excludes the axes (as in `AxisComponent`)
    const plotWidth = width - (is2d ? ySize + 10 : 0);
    const plotHeight = height - (is2d ? xSize + 10 : 0) - 2;

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
    useZoomOnFilter(data, plotWidth, plotHeight);

    const { scatterProps, selectionLayer } = useSpatialLayers();
    // this is now somewhat able to render for "2d", pending further tweaks
    //! beware unproject from here is not what we want, should review
    const { scatterplotLayer, greyScatterplotLayer, getTooltip, setScatterKeyboardActive } = scatterProps;

    const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();

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
        return config.dimension === "2d"
            ? new OrthographicView({
                  id: `scatterplot-view-${id}`,
                  // controller: true,
                  width: chartWidth,
                  height: chartHeight,
                  x: 0,
                  y: 0,
                  flipY: false,
              })
            : new OrbitView({
                  id: `scatterplot-view-${id}`,
                  // controller: true,
                  width: chartWidth,
                  height: chartHeight,
                  x: 0,
                  y: 0,
              });
    }, [chartWidth, chartHeight, config.dimension, id]);

    //! deck doesn't like it if we change the layers array - better to toggle visibility
    const layers = [
        greyScatterplotLayer,
        scatterplotLayer,
        gateDisplayLayer,
        selectionLayer,
        gateLabelLayer,
        axisLinesLayer,
    ].filter((layer): layer is NonNullable<typeof layer> => layer != null);
    const { device, failedOnWebGPU, deviceProps, onDeviceInitialized, onError } = useScatterDeckDevice(
        layers.map((layer) => layer.id),
    );
    // deck picks synchronously only on WebGL; on WebGPU the hover pick is all there is
    const pickAlternates = device.status === "ready" && device.backend === "webgl";
    const deckLayers = useLayersForDevice(layers, device);

    const deckContainerRef = useRef<HTMLDivElement | null>(null);
    const deckRef = useRef<DeckGLRef | null>(null);
    const chart = useChart<DeckScatterConfig, DeckScatterReact>();
    useEffect(() => {
        chart.getPlotCanvas = () => {
            const deck = deckRef.current?.deck;
            // drawing now lets the canvas be read before the browser clears it
            deck?.redraw("image export");
            // WebGPU work is submitted at the end of a frame, which this is not part of
            deck?.getCanvasContext()?.device?.submit();
            return deck?.getCanvas() ?? null;
        };
        return () => {
            chart.getPlotCanvas = undefined;
        };
    }, [chart]);
    const getTooltipContent = useCallback(
        (info: PickingInfo) => {
            const richInfo = getPickingInfoWithAlternates(info, pickAlternates ? deckRef.current?.deck : null);
            return getCombinedScatterTooltip(richInfo, {
                gateDisplayLayerId: gateDisplayLayer?.id,
                gateLabelLayerId: gateLabelLayer?.id,
                getPointTooltip: getTooltip,
            });
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
                    <ScatterDeckBoundary device={device} failedOnWebGPU={failedOnWebGPU}>
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
                    </ScatterDeckBoundary>
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

export default observer(function DeckScatterChart() {
    const chart = useChart<DeckScatterConfig, DeckScatterReact>();
    const [hoveredContour, setHoveredContour] = useState<HoveredContour | null>(null);
    // in order for SelectionOverlay to work, we need to review how our implementation works
    // vs useScatterplotLayer in spatial_context.tsx
    return (
        <SpatialAnnotationProvider chart={chart} hoveredContour={hoveredContour} setHoveredContour={setHoveredContour}>
            <SelectionOverlay />
            {isDensityGridShown(chart.config) ? (
                <DeckDensityGridComponent />
            ) : (
                <DeckScatter setHoveredContour={setHoveredContour} />
            )}
        </SpatialAnnotationProvider>
    );
});
