import type { Matrix4 } from "@math.gl/core";
import type { OrbitViewState, OrthographicViewState, PickingInfo } from "@deck.gl/core";
import { useChart } from "./context";
import {
    useChartID,
    useChartSize,
    useConfig,
    useFieldSpecs,
    useFilteredIndices,
    useParamColumns,
} from "./hooks";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { getVivId } from "./components/avivatorish/MDVivViewer";
import { useMetadata } from "./components/avivatorish/state";
import type { ViewState } from "./components/VivScatterComponent";
import SpatialLayer from "@/webgl/SpatialLayer";
import {
    ScatterSquareExtension,
    ScatterDensityExension,
} from "../webgl/ScatterDeckExtension";
import { useHighlightedIndex } from "./selectionHooks";
import { type DualContourLegacyConfig, useLegacyDualContour } from "./contour_state";
import type { ColumnName, FieldName } from "@/charts/charts";
import type { FeatureCollection } from "@turf/helpers";
import type { BaseConfig } from "@/charts/BaseChart";
import type { FieldSpec, FieldSpecs } from "@/lib/columnTypeHelpers";
import { getEmptyFeatureCollection } from "./deck_state";

//!!! temporary fix for tsgo preview compatibility
// import type { TooltipContent } from "@deck.gl/core/dist/lib/tooltip";
type TooltipContent = null | string | {
    text?: string;
    html?: string;
    className?: string;
    style?: Partial<CSSStyleDeclaration>;
};

export type TooltipConfig = {
    tooltip: {
        show: boolean;
        column?: FieldSpecs | FieldSpec;
    };
};
export type AxisConfig = {
    size: number;
    tickfont: number;
    rotate_labels: boolean;
};
export type AxisConfig2D = {
    x: AxisConfig;
    y: AxisConfig;
};
export type CategoryFilter = {
    column: ColumnName;
    category: string | string[];
    // consider properties like 'invert' or 'exclude', or 'color'...
};
//viewState should be persisted... maybe a way of saving different snapshots?
//now using MobX for this, in a way that's common to viv & non-viv charts
export type ScatterPlotConfig = {
    course_radius: number;
    radius: number;
    opacity: number;
    // color_by: ColumnName;
    // color_legend: {
    //     display: boolean;
    //     // todo: add more options here...
    // };
    category_filters: Array<CategoryFilter>;
    on_filter: "hide" | "grey", //todo...
    zoom_on_filter: boolean;
    point_shape: "circle" | "square" | "gaussian";
    dimension: "2d" | "3d";
    selectionFeatureCollection: FeatureCollection;
} & TooltipConfig & DualContourLegacyConfig & BaseConfig;

export type ScatterPlotConfig2D = ScatterPlotConfig & {
    dimension: "2d";
    axis: AxisConfig2D;
    viewState: OrthographicViewState;
};
export type ScatterPlotConfig3D = ScatterPlotConfig & {
    dimension: "3d";
    viewState: OrbitViewState;
};

export const scatterDefaults: Omit<ScatterPlotConfig, "id" | "legend" | "size" | "title" | "type" | "param"> = {
    //default 10 may be closer to original charts, not necessarily the best
    course_radius: 1,
    radius: 10,
    opacity: 1,
    // color_by: null,
    // color_legend: {
    //     display: false,
    // },
    tooltip: {
        show: false,
    },
    category_filters: [],
    zoom_on_filter: false,
    point_shape: "circle",
    contour_fill: false,
    contour_bandwidth: 0.1,
    contour_intensity: 1,
    contour_opacity: 0.5,
    contour_fillThreshold: 2,
    dimension: "2d",
    on_filter: "hide", //safer in case of large datasets
    // todo omit this so we can have better HMR...
    selectionFeatureCollection: getEmptyFeatureCollection(),
};

export const scatterAxisDefaults: AxisConfig2D = {
    x: {
        rotate_labels: false,
        size: 20,
        tickfont: 10,
    },
    y: {
        rotate_labels: false,
        size: 40,
        tickfont: 10,
    },
}

export function useRegionScale() {
    const metadata = useMetadata();
    const chart = useChart();
    const regionScale = chart.dataStore.regions?.scale;
    const regionUnit = chart.dataStore.regions?.scale_unit;

    //see also getPhysicalScalingMatrix
    //- consider state, matrices for image, scatterplot/other layers, and options to manipulate them
    //MDVProject.set_region_scale assumes that all regions have the same scale?
    if (!metadata) return 1 / regionScale; // might want to start using this with non-image data that has real units
    if (!("Pixels" in metadata)) return 1/regionScale;
    const { Pixels } = metadata;
    if (!Pixels.PhysicalSizeX) return 1 / regionScale;
    if (Pixels.PhysicalSizeXUnit !== regionUnit)
        console.warn(
            `physical size unit mismatch ${Pixels.PhysicalSizeXUnit} !== ${regionUnit}`,
        );
    // if (!Pixels.PhysicalSizeX) throw new Error("missing physical size");
    const scale = Pixels.PhysicalSizeX / regionScale;
    return Number.isFinite(scale) ? scale : 1;
}

/**
 * This hook is used to fit the scatterplot to the data when data filter changes.
 * 
 * It can be a bit janky when reacting to changes originating from the same view,
 * we should consider a better approach.
 */
function useZoomOnFilter(modelMatrix: Matrix4) {
    const config = useConfig<ScatterPlotConfig>();
    const data = useFilteredIndices();
    const [cx, cy] = useParamColumns();
    const [chartWidth, chartHeight] = useChartSize(); //not sure we want this, potentially re-rendering too often...
    // not using as dependency for scaling viewState to data - we don't want to zoom as chart size changes
    // (at least for now - may consider making this configurable / testing it out)

    const [viewState, setViewState] = useState<ViewState | null>(); //we should consider how this interacts with viv ViewerStore
    useEffect(() => {
        if (!config.zoom_on_filter) return;
        if (data.length === 0) {
            setViewState(null);
            return;
        }
        // Step 1: Calculate the bounding box of the data
        // with a loop to avoid spread operator & stack overflow etc
        let minX = Number.POSITIVE_INFINITY;
        let maxX = Number.NEGATIVE_INFINITY;
        let minY = Number.POSITIVE_INFINITY;
        let maxY = Number.NEGATIVE_INFINITY;
        for (let i = 0; i < data.length; i++) {
            const d = data[i];
            try {
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
            } catch (e) {
                console.error("error calculating bounding box", e);
                return;
            }
        }

        // Step 2: Calculate the center of the bounding box
        // - take into account transform matrices
        // could review using scatterplotLayer.project / projectPosition
        // quicker to do once than for every point
        [minX, minY] = modelMatrix.transformAsPoint([minX, minY, 0]);
        [maxX, maxY] = modelMatrix.transformAsPoint([maxX, maxY, 0]);

        const centerX = (minX + maxX) / 2;
        const centerY = (minY + maxY) / 2;

        // modelMatrix.transformPoint([centerX, centerY, 0]);

        // Step 3: Calculate the zoom level
        // see viv's getDefaultInitialViewState for reference, especially when it comes to 3D in the future
        // should allow some padding around the edges, based on radius, and/or some config value
        const { min, log2 } = Math;
        const maxZoom = 5;
        const trueSelectionWidth = maxX - minX; // * scale;
        const trueSelectionHeight = maxY - minY; // * scale;
        const xZoom = log2(chartWidth / trueSelectionWidth);
        const yZoom = log2(chartHeight / trueSelectionHeight);
        const zoomBackOff = 0.05;
        const zoom = min(maxZoom, min(xZoom, yZoom) - zoomBackOff);
        // Step 4: Set the view state that will be picked up by the user of this hook
        // may want to use viv viewerStore - but we don't want scatterplot to depend on that
        setViewState({
            target: [centerX, centerY, 0],
            zoom: zoom,
            // minZoom: -10,
            maxZoom,
            transitionDuration: 400, // Smooth transition
            transitionEasing: (x: number) => -(Math.cos(Math.PI * x) - 1) / 2, //https://easings.net/#easeInOutSine
            // transitionInterpolator: new FlyToInterpolator({speed: 1}), //applicable for MapState - latitude is required
        });
    }, [
        data,
        cx,
        cy,
        chartHeight,
        chartWidth,
        config.zoom_on_filter,
        modelMatrix.transformAsPoint,
    ]);
    return viewState;
}

/**
 * If we are in a 'region' then we might have some notion of physical size.
 * If we are rendering abstract data, then we should probably have a normalised size
 * with respect to the domain of the data.
 * If we have very different scales in x and y, then we probably want to also scale each axis
 * rather than assuming a 1:1 ratio - somewhat related concern... also need to consider that
 * we still want circular points to be circular...
 *
 * Also note that we would like in future to be able to set the size of individual points
 * based on some property of the data - but for now we just return a number.
 */
export function useScatterRadius() {
    const config = useConfig<ScatterPlotConfig>();
    const params = useParamColumns();
    const [cx, cy] = params;
    const scale = useRegionScale();
    const { radius, course_radius } = config;
    //todo more clarity on radius units - but large radius was causing big problems after deck upgrade
    // this is reasonably ok looking, but even for abstract data it should really relate to axis labels
    // (which implies that if we have a warped aspect ratio but making circles circular, they will be based on one or other axis)
    // see DensityPlot.js for an example of how this has been done there:
    // y_scale = [0, 400 / whRatio];
    // x_scale = [0, 400];

    const safeScale = scale || 1; // avoid /0 (although - I don't _think_ useRegionScale() would return 0).
    const radiusScale = (radius * course_radius) / safeScale;
    return useMemo(() => {
        if (cx.minMax && cy.minMax) {
            const xRange = cx.minMax[1] - cx.minMax[0];
            const yRange = cy.minMax[1] - cy.minMax[0];
            const r = 10000 / Math.max(xRange, yRange);
            return radiusScale / r;
        }
        return radiusScale;
    }, [cx.minMax, cy.minMax, radiusScale, cx, cy]);
}

// type Tooltip = (PickingInfo) => string;
export type P = [number, number];
/**
 * ! in its current form, this hook is only called by `useCreateSpatialAnnotationState`
 * in future we may want to be able to have different arrangement of layers & rework this.
 *
 * As of now, charts with appropriate spatial context can call `useSpatialLayers()` at any point
 * to access the scatterplot layer, and the tooltip function.
 */
export function useScatterplotLayer(modelMatrix: Matrix4, hoveredFieldId?: FieldName | null) {
    const id = useChartID();
    const chart = useChart();
    const colorBy = (chart as any).colorBy;
    const config = useConfig<ScatterPlotConfig>();

    const { opacity } = config;
    const radiusScale = useScatterRadius();

    const data = useFilteredIndices();
    //! not keen on third param potentially being either contourParameter or cz
    // n.b. Viv version already has config.contourParameter (maybe should be densityParameter)
    // param[2] is set to the same value for some kind of backward compatibility?
    // or as the result of still using old BaseChart.init()
    // const [cx, cy, contourParameter] = useParamColumns();
    const params = useParamColumns();
    const [cx, cy, cz] = params;
    const scale = useRegionScale();
    const hoverInfoRef = useRef<PickingInfo | null>(null);
    const highlightedIndex = useHighlightedIndex();
    // const [highlightedObjectIndex, setHighlightedObjectIndex] = useState(-1);
    const getLineWidth = useCallback(
        (i: unknown) => {
            if (typeof i !== "number") throw new Error("expected index");
            return i === highlightedIndex ? (0.2 * radiusScale) / scale : 0.0;
        },
        [radiusScale, highlightedIndex, scale],
    );
    const contourLayers = useLegacyDualContour(hoveredFieldId);

    // todo - Tooltip should be a separate component
    // would rather not even need to call a hook here, but just have some
    // `state.tooltip.column` which would have a column object...
    // but this isn't really all that bad, so maybe we can stick with it.
    const tooltipCols = useFieldSpecs(config.tooltip.column);
    const getTooltipVal = useCallback(
        (i: number) => {
            // if (!tooltipCol?.data) return '#'+i;
            if (!tooltipCols) return null;
            // return tooltipCols.getValue(data[i]);
            return tooltipCols.map((col) => {
                    return `<strong>${col.name}:</strong> ${col.data ? col.getValue(data[i]) : "loading..."}`;
            });
        },
        [tooltipCols, data],
    );
    const getTooltip = useCallback(
        //todo nicer tooltip interface (and review how this hook works)
        () => {
            if (!config.tooltip.show) return null;
            if (!config.tooltip.column) return null;
            // testing reading object properties --- pending further development (for GeoJSON layer in particular)
            // also consider some other things like transcripts / stats heatmap etc...
            // (not hardcoding DN property etc)
            // if (object && object?.properties?.DN) return `DN: ${object.properties.DN}`;
            const hoverInfo = hoverInfoRef.current;
            if (!hoverInfo || hoverInfo.index === -1) return null;
            const tooltipVal = getTooltipVal(hoverInfo.index);
            if (!tooltipVal) return null;
            const tooltip: TooltipContent = {
                //todo - this should be in a popper / should follow useOuterConainer...
                //also should understand if mouse has left deck.gl canvas & hide tooltip
                //or maybe we actually use something else for rendering the tooltip?
                html: `<div>${tooltipVal.join("<br/>")}</div>`,
            }
            return tooltip;
        },
        [getTooltipVal, config.tooltip.show, config.tooltip.column],
    );

    // const { modelMatrix, setModelMatrix } = useScatterModelMatrix();
    const viewState = useZoomOnFilter(modelMatrix);
    const { point_shape } = config;

    // could probably bring this more into SpatialLayer...
    const extensions = useMemo(() => {
        if (point_shape === "circle") return [];
        if (point_shape === "gaussian") return [new ScatterDensityExension()];
        return [new ScatterSquareExtension()];
    }, [point_shape]);
    const scatterplotLayer = useMemo(() => {
        const is3d = config.dimension === "3d";
        return new SpatialLayer({
            //new
            // loaders //<< this will be interesting to learn about
            id: `scatter_${getVivId(`${id}detail-react`)}`, // should satisfy VivViewer, could make this tidier
            data,
            opacity,
            radiusScale,
            billboard: true,
            // antialiasing has artefacts on edges in 3d.
            // may also consider option to turn off antialiasing for performance?
            // (if it even has a significant impact)
            antialiasing: !is3d,
            getFillColor: colorBy ?? [55, 126, 184],
            // if we want different radii for different points, this is the place
            // ^^ except that we probably want each density layer to have its own radius accessor
            // getRadius: 1 / scale,
            // todo review buffer data / accessors / filters...
            getPosition: (i: unknown, { target }) => {
                if (typeof i !== "number") throw new Error("expected index");
                target[0] = cx.data[i];
                target[1] = cy.data[i];
                // this `cz?.minMax` doesn't account for the potential of cz existing,
                // but not being what we anticipated
                // this shouldn't happen now - config.param should only be coordinates.
                target[2] = cz?.minMax ? cz.data[i] : 0;
                return target as unknown as Float32Array; // deck.gl types are wrong AFAICT
            },
            modelMatrix,
            updateTriggers: {
                getFillColor: colorBy, //this is working; removing it breaks the color change...
                // modelMatrix: modelMatrix, // this is not necessary, manipulating the matrix works anyway
                // getLineWith: clickIndex, // this does not work, seems to need something like a function
                getLineWidth,
                getPosition: [cx, cy, cz],
                // getRadius: [radiusScale, scale],
                //as of now, the SpatialLayer implemetation needs to figure this out for each sub-layer.
                // getContourWeight1: config.category1,
            },
            pickable: true,
            onHover: (info) => {
                hoverInfoRef.current = info;
            },
            stroked: data.length < 1000, //todo make this configurable, and fix issue...
            // todo figure out why lineWidth 0 still shows up, particularly when zoomed out
            // can we make it have zero opacity? Seems like lineColor is rgb, not rgba...
            // >>> may need a layer extension to do this properly; may want that anyway for other reasons <<<
            getLineWidth,
            //trying to set line color to same as fill, but it makes things very muddy when zoomed out
            //getLineColor: i => i === clickIndexRef.current ? [255, 255, 255] : colorBy ?? [200, 200, 200],
            // lineColorBy...
            getLineColor: [255, 255, 255],
            // highlightedObjectIndex, // has some undesirable effects, but could be useful when better controlled
            onClick: ({ index }) => {
                // setHighlightedObjectIndex(index);
                //todo properly synchronise state with data store, allow deselection
                chart.dataStore.dataHighlighted([data[index]], chart);
                // timeout allowed us to highlight & redraw this chart before heavy blocking filter operations...
                // but now we get highlight from useHighlightedIndex() we'd need more logic to short-circuit that.
                // Really want to make the filtering async etc.
                // setTimeout(()=> chart.dataStore.dataHighlighted([data[index]], chart), 5);
            },
            transitions: {
                // this leads to weird behaviour when filter changes, looks ok when changing colorBy
                // getFillColor: {
                //     duration: 300,
                // },
            },
            // ...config, //make sure contour properties are passed through
            contourLayers,
            extensions,
        });
    }, [
        id,
        data,
        opacity,
        radiusScale,
        colorBy,
        cx,
        cy,
        cz,
        modelMatrix,
        extensions,
        chart,
        getLineWidth,
        contourLayers,
        config.dimension,
    ]);
    // this should take into account axis margins... not chart.contentDiv,
    // but the actual area where the scatterplot is rendered
    // (which may in future be a smaller region within the deck.gl canvas itself)
    const boundingClientRect = useMemo(() => {
        return chart.contentDiv.getBoundingClientRect();
    }, [chart.contentDiv.getBoundingClientRect]);
    // maybe we should rename this `unprojectMouse` or something
    const unproject = useCallback(
        (e: MouseEvent | React.MouseEvent | P) => {
            // never liked this... and then it was actually causing an infinite loop
            // if (!currentLayerHasRendered || !scatterplotLayer.internalState)
            //     throw new Error("scatterplotLayer not ready");
            // also not massively keen on try/catch, but much better the dodgy currentLayerHasRendered
            try {
                if (Array.isArray(e)) {
                    e = { clientX: e[0], clientY: e[1] } as MouseEvent;
                }
                const r = boundingClientRect;
                const x = e.clientX - r.left;
                const y = e.clientY - r.top;
                const p = scatterplotLayer.unproject([x, y]) as P;
                //still need to reason better about transforms...
                // const scale = 1; //this was only right when Pixels.PhysicalSizeX === regions.scale
                // const p2 = p.map(v => v * scale) as P;
                const m = modelMatrix.invert();
                const p3 = m.transform(p) as P;
                m.invert();
                return p3;
            } catch (e) {
                //console.error("unproject error", e);
                return [0, 0];
            }
        },
        [
            scatterplotLayer,
            modelMatrix,
            boundingClientRect,
        ],
    );
    return useMemo(
        () => ({
            scatterplotLayer,
            getTooltip,
            modelMatrix,
            viewState,
            unproject,
        }),
        [
            scatterplotLayer,
            getTooltip,
            modelMatrix,
            viewState,
            unproject,
        ],
    );
}
