import { type FilterContext, OrbitView, OrthographicView, type PickingInfo } from "@deck.gl/core";
import { LineLayer } from "@deck.gl/layers";
import DeckGL, { type DeckGLRef } from "@deck.gl/react";
import { action, toJS } from "mobx";
import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useId, useMemo, useRef, useState } from "react";
import type { LoadedDataColumn } from "@/charts/charts";
import type { MonkeyPatchEditableGeoJsonLayer } from "@/lib/deckMonkeypatch";
import { getPickingInfoWithAlternates } from "@/lib/deckPicking";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import { GatedOrbitController, GatedOrthographicController } from "@/lib/viewGestureGate";
import { useChart } from "../context";
import {
    type DualContourLegacyConfig,
    type HoveredContour,
    useCategoryContourLegend,
    useFieldContourLegend,
} from "../contour_state";
import { getEmptyFeatureCollection } from "../deck_state";
import { useChartSize, useConfig, useFilteredIndices, useParamColumns } from "../hooks";
import { useDensityGridCells } from "../hooks/useDensityGridCells";
import useGateLayers from "../hooks/useGateLayers";
import { useOuterContainerDeckTooltip } from "../hooks/useOuterContainerDeckTooltip";
import { useRebindDeckEvents } from "../hooks/useRebindDeckEvents";
import { isSameOrbitView, ORBIT_FOVY, type OrbitSelectionView, orbitSelectionViewOf } from "../orbitView";
import { SpatialAnnotationProvider, useSpatialLayers } from "../spatial_context";
import AxisComponent from "./AxisComponent";
import DeckDensityGridComponent from "./DeckDensityGridComponent";
import type { DeckScatterConfig, DeckScatterReact } from "./DeckScatterReactWrapper";
import { isDensityGridShown } from "./densityGridUtils";
import { useContourLegend } from "./FieldContourLegend";
import SelectionOverlay from "./SelectionOverlay";
import { ScatterDeckBoundary, useLayersForDevice, useScatterDeckDevice, useZoomOnFilter } from "./scatterDeckShared";
import ViewGestureHint from "./ViewGestureHint";

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
const EMPTY_SELECTION = getEmptyFeatureCollection();

const selectionViewId = (id: string) => `selection-screen-${id}`;

/**
 * In 3D, selections are drawn on the screen rather than in the data: a view in pixels over the orbit
 * view holds the selection layer, and nothing else. A selection is what was inside its outline from
 * the camera it was drawn with, so its outline is only shown while the camera is still there; it
 * stays applied (and saved) when the camera moves, until it is cleared or another is drawn.
 */
function useScreenSelection(
    id: string,
    selectionLayer: MonkeyPatchEditableGeoJsonLayer,
    width: number,
    height: number,
) {
    const chart = useChart<DeckScatterConfig, DeckScatterReact>();
    const config = useConfig<DeckScatterConfig>();
    const is3d = config.dimension === "3d";
    const orbitViewState = config.dimension === "3d" ? config.viewState : null;
    const orbitAxis = config.dimension === "3d" ? (config.orbit_axis ?? "Z") : "Z";
    const current = useMemo<OrbitSelectionView | null>(() => {
        if (!orbitViewState) return null;
        const { target, zoom, rotationOrbit = 0, rotationX = 0 } = orbitViewState;
        return {
            target: [target[0], target[1], target[2] ?? 0],
            zoom,
            rotationOrbit,
            rotationX,
            width,
            height,
            orbitAxis,
        };
    }, [orbitViewState, width, height, orbitAxis]);
    useEffect(() => {
        chart.getOrbitSelectionView = () => current;
        return () => {
            chart.getOrbitSelectionView = undefined;
        };
    }, [chart, current]);

    const feature = config.selectionFeatureCollection.features[0];
    const drawnIn = useMemo(() => orbitSelectionViewOf(toJS(feature?.properties)), [feature]);
    const shown = !is3d || (current && drawnIn && isSameOrbitView(drawnIn, current));
    const layer = useMemo(
        () => (shown ? selectionLayer : selectionLayer.clone({ data: EMPTY_SELECTION, selectedFeatureIndexes: [] })),
        [shown, selectionLayer],
    );
    const screenViewId = selectionViewId(id);
    const layerFilter = useCallback(
        ({ layer: l, viewport }: FilterContext) =>
            viewport.id === screenViewId ? l.id === selectionLayer.id : l.id !== selectionLayer.id,
        [screenViewId, selectionLayer.id],
    );
    return {
        selectionLayer: layer,
        layerFilter: is3d ? layerFilter : undefined,
        viewState: orbitViewState
            ? {
                  [`scatterplot-view-${id}`]: orbitViewState,
                  // world units are pixels from the top left
                  [screenViewId]: { target: [width / 2, height / 2, 0], zoom: 0 },
              }
            : undefined,
    };
}

const DeckScatter = observer(function DeckScatterComponent({
    setHoveredContour,
}: {
    setHoveredContour: (contour: HoveredContour | null) => void;
}) {
    const id = useId();
    const [width, height] = useChartSize();
    const [cx, cy] = useParamColumns() as LoadedDataColumn<"double">[];
    const cz = useParamColumns()[2] as LoadedDataColumn<"double">;
    const data = useFilteredIndices(); //changed to fallback to simplerFilteredIndices when filterColumn is not set
    const config = useConfig<DeckScatterConfig>();
    const contourConfig = useConfig<DualContourLegacyConfig>();
    const { viewState, dimension } = config;
    const is2d = dimension === "2d";
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

    const { scatterProps, selectionLayer, selectionProps } = useSpatialLayers();
    // this is now somewhat able to render for "2d", pending further tweaks
    //! beware unproject from here is not what we want, should review
    const { scatterplotLayer, greyScatterplotLayer, unselectedScatterplotLayer, getTooltip, setScatterKeyboardActive } =
        scatterProps;

    const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();

    const legendFields = useFieldContourLegend(contourConfig.densityFields);
    const legendCategories = useCategoryContourLegend(contourConfig);
    // drawn in the chart's legend wrapper: draggable, can be minimised or put above the plot
    useContourLegend(contourConfig.field_legend, {
        fields: legendFields,
        categories: legendCategories.items,
        categoryLabel: legendCategories.label,
        onContourHover: setHoveredContour,
    });

    const { modelMatrix } = selectionProps;
    // as the classic 3D scatter drew them: each axis across its column's range, through the middle of the data
    const axisLinesLayer = useMemo(() => {
        if (is2d || !cz) return null;
        const ranges = [cx.minMax, cy.minMax, cz.minMax];
        const mid = ranges.map(([lo, hi]) => (lo + hi) / 2);
        const axisLine = (axis: number) => ({
            sourcePosition: mid.map((m, i) => (i === axis ? ranges[i][0] : m)),
            targetPosition: mid.map((m, i) => (i === axis ? ranges[i][1] : m)),
            color: [0, 1, 2].map((i) => (i === axis ? 255 : 0)),
        });
        return new LineLayer({
            id: `axis-lines-${id}`,
            data: [axisLine(0), axisLine(1), axisLine(2)],
            getSourcePosition: (d: any) => d.sourcePosition,
            getTargetPosition: (d: any) => d.targetPosition,
            getColor: (d: any) => d.color,
            getWidth: 1,
            modelMatrix,
        });
    }, [is2d, id, cx.minMax, cy.minMax, cz?.minMax, cz, modelMatrix]);

    // we need an OrthographicView to prevent wrapping etc...
    // if in future we have subgraphs sharing a canvas, we will need to
    // make sure that the view is set up correctly for each subgraph.
    const orbitAxis = config.dimension === "3d" ? (config.orbit_axis ?? "Z") : "Z";
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
            : [
                  // first, so deck gives it the controller
                  new OrbitView({
                      id: `scatterplot-view-${id}`,
                      width: chartWidth,
                      height: chartHeight,
                      x: 0,
                      y: 0,
                      orbitAxis,
                      fovy: ORBIT_FOVY,
                  }),
                  // over it, in pixels, for drawing selections on (see `useScreenSelection`)
                  new OrthographicView({
                      id: selectionViewId(id),
                      width: chartWidth,
                      height: chartHeight,
                      x: 0,
                      y: 0,
                      controller: false,
                  }),
              ];
    }, [chartWidth, chartHeight, config.dimension, id, orbitAxis]);
    const screenSelection = useScreenSelection(id, selectionLayer, chartWidth, chartHeight);

    //! deck doesn't like it if we change the layers array - better to toggle visibility
    // gates are x / y outlines, which don't say what they enclose in 3D
    const layers = [
        greyScatterplotLayer,
        // drawn first and leaving no depth behind, so what a 3D selection caught shows through the rest
        unselectedScatterplotLayer,
        scatterplotLayer,
        is2d ? gateDisplayLayer : null,
        screenSelection.selectionLayer,
        is2d ? gateLabelLayer : null,
        axisLinesLayer,
    ].filter((layer): layer is NonNullable<typeof layer> => layer != null);
    const { device, failedOnWebGPU, deviceProps, onDeviceInitialized, onError } = useScatterDeckDevice(
        layers.map((layer) => layer.id),
    );
    // deck picks synchronously only on WebGL; on WebGPU the hover pick is all there is
    const pickAlternates = device.status === "ready" && device.backend === "webgl";
    const deckLayers = useLayersForDevice(layers, device);

    const deckContainerRef = useRef<HTMLElement | null>(null);
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
                <section
                    ref={deckContainerRef}
                    aria-label="Scatter plot"
                    style={{ width: "100%", height: "100%", outline: "none", position: "relative" }}
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
                                type: is2d ? GatedOrthographicController : GatedOrbitController,
                                // double-click finishes polygons / removes vertices in selection tools
                                doubleClickZoom: false,
                                ...controllerOptions,
                            }}
                            viewState={screenSelection.viewState ?? viewState}
                            // initialViewState={viewState} //consider not using react state for this
                            views={view}
                            layerFilter={screenSelection.layerFilter}
                            onViewStateChange={(v) => {
                                // the selection view does not move
                                if (v.viewId === selectionViewId(id)) return;
                                action(() => (config.viewState = v.viewState))();
                            }}
                            getTooltip={getPortalTooltip}
                            getCursor={({ isDragging }) => {
                                return isDragging ? "grabbing" : "crosshair";
                            }}
                        />
                    </ScatterDeckBoundary>
                </section>
            </AxisComponent>
            <ViewGestureHint />
            {tooltipPortal}
        </>
    );
});

export default observer(function DeckScatterChart() {
    const chart = useChart<DeckScatterConfig, DeckScatterReact>();
    const [hoveredContour, setHoveredContour] = useState<HoveredContour | null>(null);
    const grid = useDensityGridCells();
    // in order for SelectionOverlay to work, we need to review how our implementation works
    // vs useScatterplotLayer in spatial_context.tsx
    return (
        <SpatialAnnotationProvider chart={chart} hoveredContour={hoveredContour} setHoveredContour={setHoveredContour}>
            <SelectionOverlay />
            {/* with no field to make a cell of (none chosen, nothing highlighted for an active link,
                or still loading) the grid would be empty, so show the one plot, without an overlay */}
            {isDensityGridShown(chart.config) && grid.cells.length > 0 ? (
                <DeckDensityGridComponent {...grid} />
            ) : (
                <DeckScatter setHoveredContour={setHoveredContour} />
            )}
        </SpatialAnnotationProvider>
    );
});
