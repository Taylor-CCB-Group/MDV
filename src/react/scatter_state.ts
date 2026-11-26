import type { LayerExtension, OrbitViewState, OrthographicViewState, PickingInfo } from "@deck.gl/core";
import type { SimpleFeatureCollection } from "@deck.gl-community/editable-layers";
import type { Matrix4 } from "@math.gl/core";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { BaseConfig } from "@/charts/BaseChart";
import type { ColumnName, DataColumn, DataType } from "@/charts/charts";
import type { FieldSpec, FieldSpecs } from "@/lib/columnTypeHelpers";
import { isDatatypeCategorical } from "@/lib/utils";
import { escapeHtml } from "@/utilities/Utilities";
import { MdvDataFilterExtension } from "@/webgl/MdvDataFilterExtension";
import { MdvScatterplotLayer } from "@/webgl/MdvScatterplotLayer";
import SpatialLayer, { DensityOverlayLayer } from "@/webgl/SpatialLayer";
import { ScatterDensityExension, ScatterSquareExtension } from "../webgl/ScatterDeckExtension";
import { type AxisScaling, radiusReferenceRange, radiusReferenceRange3D } from "./axisScaling";
import { getVivId } from "./components/avivatorish/MDVivViewer";
import { useMetadata } from "./components/avivatorish/state";
import type { ViewState } from "./components/VivScatterComponent";
import { useChart } from "./context";
import {
    type ContourLayerProps,
    type DualContourLegacyConfig,
    type HoveredContour,
    useLegacyDualContour,
} from "./contour_state";
import { getEmptyFeatureCollection } from "./deck_state";
import {
    type FilterOwner,
    useChartID,
    useChartSize,
    useConfig,
    useFieldSpec,
    useFieldSpecs,
    useOwnedFilteredIndices,
    useParamColumns,
} from "./hooks";
import type { AxisScales, OrbitAxis } from "./orbitView";
import { useHighlightedIndices, useHighlightRows } from "./selectionHooks";

//!!! temporary fix for tsgo preview compatibility
// import type { TooltipContent } from "@deck.gl/core/dist/lib/tooltip";
type TooltipContent =
    | null
    | string
    | {
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
/**
 * Hides rows outside `category` of `column` from this chart (`null` when unset).
 * Read by `useChartScopeFilterPredicate`; it does not filter other charts.
 */
export type BackgroundFilter = {
    column: ColumnName | null;
    category: string | string[];
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
    background_filter?: BackgroundFilter | null;
    on_filter: "hide" | "grey"; //todo...
    zoom_on_filter: boolean;
    /** Show every selection tool; when false the toolbar shows one button that switches between Pan and the last drawing tool. */
    selection_toolbar_open: boolean;
    /** 2D only: mark the pointer with a crosshair; absent is on. null only keeps the key observable. */
    crosshair?: boolean | null;
    point_shape: "circle" | "square" | "gaussian";
    /** 2D only; absent is `equal`. See `AxisScaling`. */
    axis_scaling?: AxisScaling;
    dimension: "2d" | "3d";
    selectionFeatureCollection: SimpleFeatureCollection;
} & TooltipConfig &
    DualContourLegacyConfig &
    BaseConfig;

export type ScatterPlotConfig2D = ScatterPlotConfig & {
    dimension: "2d";
    axis: AxisConfig2D;
    viewState: OrthographicViewState;
};
export type ScatterPlotConfig3D = ScatterPlotConfig & {
    dimension: "3d";
    viewState: OrbitViewState;
    /** The data axis that points up and that the view orbits round; absent is Z. */
    orbit_axis?: OrbitAxis;
    /** How much to stretch the x, y and z axes, as the classic 3D scatter did; absent is [1, 1, 1]. */
    axis_scales?: AxisScales;
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
    selection_toolbar_open: true,
    point_shape: "circle",
    contourParameter: undefined,
    category1: [],
    category2: [],
    contour_fill: false,
    contour_bandwidth: 0.1,
    contour_intensity: 1,
    contour_opacity: 0.5,
    contour_fillThreshold: 2,
    densityFields: [],
    dimension: "2d",
    on_filter: "grey",
    // todo omit this so we can have better HMR...
    selectionFeatureCollection: getEmptyFeatureCollection(),
    field_legend: {
        display: true,
    },
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
};

export function useRegionScale() {
    const metadata = useMetadata();
    const chart = useChart();
    const regionScale = chart.dataStore.regions?.scale;
    const regionUnit = chart.dataStore.regions?.scale_unit;
    //! Fix for the scatterplot not showing suggested by cursor
    // guard against missing or invalid scale to avoid NaNs downstream
    if (!regionScale || !Number.isFinite(regionScale) || regionScale === 0) {
        return 1;
    }

    // Fix for the scatterplot not showing suggested by cursor
    // guard against missing or invalid scale to avoid NaNs downstream
    if (!regionScale || !Number.isFinite(regionScale) || regionScale === 0) {
        return 1;
    }

    //see also getPhysicalScalingMatrix
    //- consider state, matrices for image, scatterplot/other layers, and options to manipulate them
    //MDVProject.set_region_scale assumes that all regions have the same scale?
    if (!metadata) return 1 / regionScale; // might want to start using this with non-image data that has real units
    if (!("Pixels" in metadata)) return 1 / regionScale;
    const { Pixels } = metadata;
    if (!Pixels.PhysicalSizeX) return 1 / regionScale;
    if (Pixels.PhysicalSizeXUnit !== regionUnit)
        console.warn(`physical size unit mismatch ${Pixels.PhysicalSizeXUnit} !== ${regionUnit}`);
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
function useZoomOnFilter(modelMatrix: Matrix4, data: Uint32Array) {
    const config = useConfig<ScatterPlotConfig>();
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
    }, [data, cx, cy, chartHeight, chartWidth, config.zoom_on_filter, modelMatrix.transformAsPoint]);
    return viewState;
}

/**
 * A radius, so 0.5 keeps the smallest point about 1px across; antialiasing fades it below that.
 * Suspect we could get away with reducing it further - on dense plots subpixel sizes can look much better.
 * This is a conservatively high value for the time being.
 */
export const SCATTER_RADIUS_MIN_PIXELS = 0.5;

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
 *
 * The result is in data units, ~`radius * course_radius / 10000` of the data range, so on screen
 * it scales with chart size and zoom: in a small chart even the largest slider value can be
 * sub-pixel. Layers using it should also set `radiusMinPixels: SCATTER_RADIUS_MIN_PIXELS`.
 */
export function useScatterRadius() {
    const { radius, course_radius, axis_scaling, dimension } = useConfig<ScatterPlotConfig>();
    const [width, height] = useChartSize();
    const scaling = dimension === "2d" ? (axis_scaling ?? "equal") : "equal";
    const params = useParamColumns();
    const [cx, cy, cz] = params;
    const scale = useRegionScale();
    //todo more clarity on radius units - but large radius was causing big problems after deck upgrade
    // this is reasonably ok looking, but even for abstract data it should really relate to axis labels
    // (which implies that if we have a warped aspect ratio but making circles circular, they will be based on one or other axis)
    // see DensityPlot.js for an example of how this has been done there:
    // y_scale = [0, 400 / whRatio];
    // x_scale = [0, 400];

    // avoid /0 (although - I don't _think_ useRegionScale() would return 0). A 3D chart's axes need
    // not be the region's x / y, so it is not drawn at the region's scale (see `useScatterModelMatrix`).
    const safeScale = dimension === "3d" ? 1 : scale || 1;
    const radiusScale = (radius * course_radius) / safeScale;
    return useMemo(() => {
        if (cx.minMax && cy.minMax) {
            // apparently `minMax` is sometimes `maxMin` - so adding abs.
            // how does this come to be? is it a spatial conversion thing?
            const xRange = Math.abs(cx.minMax[1] - cx.minMax[0]);
            const yRange = Math.abs(cy.minMax[1] - cy.minMax[0]);
            const reference =
                dimension === "3d" && cz?.minMax
                    ? radiusReferenceRange3D([xRange, yRange, cz.minMax[1] - cz.minMax[0]])
                    : radiusReferenceRange(xRange, yRange, { width, height }, scaling);
            const r = 10000 / reference;
            return radiusScale / r;
        }
        return radiusScale;
    }, [cx.minMax, cy.minMax, cz?.minMax, radiusScale, cx, cy, width, height, scaling, dimension]);
}

/** How visible the rows a 3D selection left out are, relative to the points it caught. */
const UNSELECTED_OPACITY = 0.08;

/**
 * Which of `rows` the chart's own filter (`owner`) leaves out, as a mask over every row (1: left
 * out), or null without an active filter. A chart draws the rows its own filter leaves out (so a
 * selection can be redrawn); this tells them apart, for a chart that shows its selection by fading
 * the rest. A mask rather than two lists of rows, so the layers keep their `data` (and deck its
 * attributes) and only the filter values change.
 */
export function ownSelectionMask(rows: Uint32Array, size: number, owner?: FilterOwner | null) {
    if (!owner?.filterMethod) return null;
    const local = owner.getLocalFilter();
    const mask = new Uint8Array(size);
    let any = false;
    for (let i = 0; i < rows.length; i++) {
        const row = rows[i];
        // bit 1: left out by this filter (2 is the background filter)
        if (local[row] & 1) {
            mask[row] = 1;
            any = true;
        }
    }
    return any ? mask : null;
}

// type Tooltip = (PickingInfo) => string;
export type P = [number, number];

export function useShouldFilterNaN() {
    const chart = useChart();
    const colorBySpec = chart.config.color_by;
    const colorByField: FieldSpec | undefined = colorBySpec
        ? typeof colorBySpec === "string"
            ? colorBySpec
            : "column" in colorBySpec
              ? colorBySpec.column?.field
              : undefined
        : undefined;
    const colorColumn = useFieldSpec(colorByField);
    const hideMissing = !!chart.config.hideMissing;
    const fallbackOnZero = !!chart.config.fallbackOnZero;
    const isNumericColor = !!(colorColumn && !isDatatypeCategorical(colorColumn.datatype));
    // We only apply the deck.gl DataFilterExtension when `hideMissing` is enabled and the
    // color column is numeric. `fallbackOnZero` then extends what counts as \"missing\"
    // inside the filter (zeros as well as NaN/Inf), but does not by itself trigger filtering.
    const shouldFilter = hideMissing && isNumericColor;
    return { shouldFilter, colorColumn, fallbackOnZero };
}

export function getMissingColorFilterValue(
    index: number,
    colorColumn: DataColumn<DataType> | undefined,
    shouldFilterMissing: boolean,
    fallbackOnZero: boolean,
): 0 | 1 {
    if (!shouldFilterMissing || !colorColumn?.data) return 1;
    const value = colorColumn.data[index];
    if (!Number.isFinite(value)) return 0;
    if (fallbackOnZero && value === 0) return 0;
    return 1;
}

/** Stable, so a points layer drawn without density keeps one `contourLayers` identity. */
const NO_CONTOURS: ContourLayerProps[] = [];

/**
 * ! in its current form, this hook is only called by `useCreateSpatialAnnotationState`
 * in future we may want to be able to have different arrangement of layers & rework this.
 *
 * As of now, charts with appropriate spatial context can call `useSpatialLayers()` at any point
 * to access the scatterplot layer, and the tooltip function.
 */
export function useScatterplotLayer(
    modelMatrix: Matrix4,
    hoveredContour?: HoveredContour | null,
    filterOwner?: FilterOwner | null,
    /** Draw density contours as `densityLayer`, apart from the points, rather than inside them. */
    separateDensityLayer = false,
) {
    const id = useChartID();
    const chart = useChart();
    const colorBy = (chart as any).colorBy;
    const config = useConfig<ScatterPlotConfig>();

    const { opacity } = config;
    const radiusScale = useScatterRadius();

    const { aggregateFilteredRows, ownerVisibleRows, isExternallyFiltered } = useOwnedFilteredIndices(filterOwner);
    // a 3D selection's outline is only shown from where it was drawn, so its rows are what show it
    const is3d = config.dimension === "3d";
    const data = ownerVisibleRows;
    const params = useParamColumns();
    const rowCount = params[0].data.length;
    const leftOut = useMemo(() => {
        // `aggregateFilteredRows` changes with every filter, this chart's own included
        aggregateFilteredRows;
        return is3d ? ownSelectionMask(ownerVisibleRows, rowCount, filterOwner) : null;
    }, [ownerVisibleRows, aggregateFilteredRows, is3d, filterOwner, rowCount]);
    const [cx, cy, cz] = params;
    const scale = useRegionScale();
    const hoverInfoRef = useRef<PickingInfo | null>(null);
    const scatterKeyboardActiveRef = useRef(false);
    const spaceHighlightActiveRef = useRef(false);
    const spaceHighlightModeRef = useRef<"add" | "remove">("add");
    const lastSpaceHighlightedRowRef = useRef(-1);
    const highlightedIndices = useHighlightedIndices();
    const highlightRows = useHighlightRows();
    const getModifierState = useCallback(
        (
            event?: {
                shiftKey?: boolean;
                ctrlKey?: boolean;
                metaKey?: boolean;
                srcEvent?: unknown;
                sourceEvent?: unknown;
            } | null,
        ) => {
            const raw = event ?? {};
            const nested =
                typeof raw === "object" && raw !== null
                    ? ((raw.srcEvent ?? raw.sourceEvent) as
                          | {
                                shiftKey?: boolean;
                                ctrlKey?: boolean;
                                metaKey?: boolean;
                                srcEvent?: unknown;
                                sourceEvent?: unknown;
                            }
                          | undefined)
                    : undefined;
            const resolved = nested ?? raw;
            return {
                add: !!resolved.shiftKey,
                toggle: !!(resolved.ctrlKey || resolved.metaKey),
            };
        },
        [],
    );
    const highlightLineWidth = useMemo(() => (0.2 * radiusScale) / scale, [radiusScale, scale]);
    const highlightedData = useMemo(() => {
        if (highlightedIndices.length === 0 || data.length === 0) return [];
        const visibleRows = new Set<number>(data as Iterable<number>);
        return highlightedIndices.filter((rowIndex) => visibleRows.has(rowIndex));
    }, [data, highlightedIndices]);
    const contourLayers = useLegacyDualContour(hoveredContour);
    const updateHighlightedRows = useCallback(
        (
            rowIndex: number,
            modifiers?: {
                add?: boolean;
                remove?: boolean;
                toggle?: boolean;
            },
        ) => {
            if (rowIndex < 0) return;
            const { add = false, remove = false, toggle = false } = modifiers ?? {};
            if (toggle) {
                const nextHighlights = new Set(highlightedIndices);
                if (nextHighlights.has(rowIndex)) nextHighlights.delete(rowIndex);
                else nextHighlights.add(rowIndex);
                highlightRows(Array.from(nextHighlights));
                return;
            }
            if (remove) {
                if (!highlightedIndices.includes(rowIndex)) return;
                highlightRows(highlightedIndices.filter((index) => index !== rowIndex));
                return;
            }
            if (add) {
                if (highlightedIndices.includes(rowIndex)) return;
                highlightRows([...highlightedIndices, rowIndex]);
                return;
            }
            highlightRows([rowIndex]);
        },
        [highlightRows, highlightedIndices],
    );
    const paintHoveredRow = useCallback(() => {
        const hoverInfo = hoverInfoRef.current;
        if (!hoverInfo || hoverInfo.index === -1) return;
        const rowIndex = data[hoverInfo.index];
        if (rowIndex === undefined || rowIndex === lastSpaceHighlightedRowRef.current) return;
        lastSpaceHighlightedRowRef.current = rowIndex;
        updateHighlightedRows(rowIndex, {
            add: spaceHighlightModeRef.current === "add",
            remove: spaceHighlightModeRef.current === "remove",
        });
    }, [data, updateHighlightedRows]);

    const { shouldFilter: shouldFilterMissing, colorColumn, fallbackOnZero } = useShouldFilterNaN();

    // todo - Tooltip should be a separate component
    // would rather not even need to call a hook here, but just have some
    // `state.tooltip.column` which would have a column object...
    // but this isn't really all that bad, so maybe we can stick with it.
    const tooltipCols = useFieldSpecs(config.tooltip.column);
    const getTooltipVal = useCallback(
        (rowIndex: number) => {
            // if (!tooltipCol?.data) return '#'+i;
            if (!tooltipCols) return null;
            // return tooltipCols.getValue(data[i]);
            return tooltipCols.map((col) => {
                // Sanitise the strings before passing
                const value = col.data ? col.getValue(rowIndex) : "loading...";
                const strValue = value === undefined || value === null ? "" : String(value);
                return `<strong>${escapeHtml(col.name)}:</strong> ${escapeHtml(strValue)}`;
            });
        },
        [tooltipCols],
    );
    const getTooltipRowIndex = useCallback(
        (info?: PickingInfo | null) => {
            const pickingInfo = info === undefined ? hoverInfoRef.current : info;
            if (!pickingInfo || pickingInfo.index === -1) return undefined;
            if (typeof pickingInfo.object === "number") return pickingInfo.object;
            return data[pickingInfo.index];
        },
        [data],
    );
    const getTooltip = useCallback(
        //todo nicer tooltip interface (and review how this hook works)
        (info?: PickingInfo | null) => {
            if (!config.tooltip.show) return null;
            if (!config.tooltip.column) return null;
            // testing reading object properties --- pending further development (for GeoJSON layer in particular)
            // also consider some other things like transcripts / stats heatmap etc...
            // (not hardcoding DN property etc)
            // if (object && object?.properties?.DN) return `DN: ${object.properties.DN}`;
            const rowIndex = getTooltipRowIndex(info);
            if (rowIndex === undefined) return null;
            const tooltipVal = getTooltipVal(rowIndex);
            if (!tooltipVal) return null;
            const tooltip: TooltipContent = {
                //todo - this should be in a popper / should follow useOuterConainer...
                //also should understand if mouse has left deck.gl canvas & hide tooltip
                //or maybe we actually use something else for rendering the tooltip?
                html: `<div>${tooltipVal.join("<br/>")}</div>`,
            };
            return tooltip;
        },
        [getTooltipRowIndex, getTooltipVal, config.tooltip.show, config.tooltip.column],
    );

    // const { modelMatrix, setModelMatrix } = useScatterModelMatrix();
    const viewState = useZoomOnFilter(modelMatrix, aggregateFilteredRows);
    const { point_shape } = config;

    const extensions = useMemo(() => {
        const exts: LayerExtension[] = [];
        if (point_shape === "gaussian") exts.push(new ScatterDensityExension());
        else if (point_shape !== "circle") exts.push(new ScatterSquareExtension());
        exts.push(new MdvDataFilterExtension({ filterSize: 1 }));
        return exts;
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
            radiusMinPixels: SCATTER_RADIUS_MIN_PIXELS,
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
            ...({
                // todo - consider lower overhead version of this.
                // future work https://deck.gl/docs/developer-guide/performance#use-binary-data
                //currently useFieldSpec is typed as DataColumn|undefined,
                //if not undefined is guaranteed to be loaded but we may change how we manage lazy-loading.
                getFilterValue: (i: number) =>
                    leftOut?.[i] ? 0 : getMissingColorFilterValue(i, colorColumn, shouldFilterMissing, fallbackOnZero),
                filterRange: [0.5, 1],
            } as any),
            updateTriggers: {
                getFillColor: colorBy,
                getPosition: [cx, cy, cz],
                getFilterValue: [colorColumn, shouldFilterMissing, fallbackOnZero, leftOut],
            },
            pickable: true,
            onHover: (info) => {
                hoverInfoRef.current = info;
                if (!spaceHighlightActiveRef.current) return;
                if (info.index === -1) {
                    lastSpaceHighlightedRowRef.current = -1;
                    return;
                }
                paintHoveredRow();
            },
            stroked: false,
            // highlightedObjectIndex, // has some undesirable effects, but could be useful when better controlled
            onClick: ({ index }, event) => {
                if (index < 0) return;
                updateHighlightedRows(data[index], getModifierState(event));
            },
            transitions: {
                // this leads to weird behaviour when filter changes, looks ok when changing colorBy
                // getFillColor: {
                //     duration: 300,
                // },
            },
            // ...config, //make sure contour properties are passed through
            contourLayers: separateDensityLayer ? NO_CONTOURS : contourLayers,
            highlightedData,
            highlightLineWidth,
            highlightRadiusScale: radiusScale * 1.15,
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
        contourLayers,
        separateDensityLayer,
        highlightedData,
        highlightLineWidth,
        config.dimension,
        shouldFilterMissing,
        colorColumn,
        fallbackOnZero,
        leftOut,
        getModifierState,
        updateHighlightedRows,
        paintHoveredRow,
    ]);

    const densityLayer = useMemo(
        () =>
            separateDensityLayer && contourLayers.length > 0
                ? new DensityOverlayLayer({
                      id: `density_${getVivId(`${id}detail-react`)}`,
                      contourLayers,
                      modelMatrix,
                      // as when drawn inside the points layer, which is pickable
                      pickable: true,
                  })
                : null,
        [separateDensityLayer, contourLayers, modelMatrix, id],
    );

    const greyOnFilter = config.on_filter === "grey";
    // deck.gl regenerates every attribute when `data` changes identity (and updates hidden
    // layers too), so this must not be a fresh object each time the layer is rebuilt
    const greyData = useMemo(() => ({ length: rowCount }), [rowCount]);
    const greyScatterplotLayer = useMemo(
        () =>
            new MdvScatterplotLayer({
                id: `scatter-grey_${getVivId(`${id}detail-react`)}`,
                data: greyData,
                opacity,
                stroked: false,
                filled: true,
                radiusScale,
                radiusMinPixels: SCATTER_RADIUS_MIN_PIXELS,
                getPosition: (_: unknown, { target, index }: { target: number[]; index: number }) => {
                    target[0] = cx.data[index];
                    target[1] = cy.data[index];
                    if (cz) target[2] = cz.data[index];
                    return target as unknown as [number, number];
                },
                getFillColor: [200, 200, 200],
                getLineColor: [0, 0, 0],
                billboard: true,
                modelMatrix,
                parameters: {
                    // portable form of WebGL's `depthTest: false`, which WebGPU rejects
                    depthCompare: "always",
                    depthWriteEnabled: false,
                },
                getFilterValue: (_: unknown, { index }: { index: number }) => {
                    if (!isExternallyFiltered(index)) return 0;
                    return getMissingColorFilterValue(index, colorColumn, shouldFilterMissing, fallbackOnZero);
                },
                filterRange: [0.5, 1],
                updateTriggers: {
                    // a full pass over every row; only worth redoing while the layer is shown
                    getFilterValue: greyOnFilter
                        ? [
                              aggregateFilteredRows,
                              isExternallyFiltered,
                              colorColumn,
                              shouldFilterMissing,
                              fallbackOnZero,
                          ]
                        : false,
                    getPosition: [cx.data, cy.data, cz?.data],
                },
                extensions: [new MdvDataFilterExtension()],
                visible: greyOnFilter,
            }),
        [
            id,
            greyData,
            cx,
            cy,
            cz,
            opacity,
            radiusScale,
            modelMatrix,
            greyOnFilter,
            aggregateFilteredRows,
            isExternallyFiltered,
            colorColumn,
            shouldFilterMissing,
            fallbackOnZero,
        ],
    );

    // the rows a 3D selection left out, faint grey; drawn first, and leaving no depth behind, so the
    // points the selection caught show through them
    const unselectedScatterplotLayer = useMemo(
        () =>
            new MdvScatterplotLayer({
                id: `scatter-unselected_${getVivId(`${id}detail-react`)}`,
                // every row, so a new selection changes only the filter values
                data: greyData,
                opacity: opacity * UNSELECTED_OPACITY,
                stroked: false,
                filled: true,
                radiusScale,
                radiusMinPixels: SCATTER_RADIUS_MIN_PIXELS,
                getPosition: (_: unknown, { target, index }: { target: number[]; index: number }) => {
                    target[0] = cx.data[index];
                    target[1] = cy.data[index];
                    target[2] = cz?.minMax ? cz.data[index] : 0;
                    return target as unknown as [number, number];
                },
                // grey, as for "Gray Out Points": colour is what the selection caught
                getFillColor: [200, 200, 200],
                billboard: true,
                modelMatrix,
                parameters: { depthWriteEnabled: false },
                getFilterValue: (_: unknown, { index }: { index: number }) =>
                    leftOut?.[index]
                        ? getMissingColorFilterValue(index, colorColumn, shouldFilterMissing, fallbackOnZero)
                        : 0,
                filterRange: [0.5, 1],
                updateTriggers: {
                    getPosition: [cx, cy, cz],
                    // a full pass over every row; only worth redoing while the layer is shown
                    getFilterValue: leftOut ? [colorColumn, shouldFilterMissing, fallbackOnZero, leftOut] : false,
                },
                extensions: [new MdvDataFilterExtension()],
                // deck prefers a layer hidden to one removed
                visible: leftOut !== null,
            }),
        [
            id,
            greyData,
            leftOut,
            opacity,
            radiusScale,
            cx,
            cy,
            cz,
            modelMatrix,
            colorColumn,
            shouldFilterMissing,
            fallbackOnZero,
        ],
    );

    const onScatterGlobalKeyDown = useCallback(
        (event: KeyboardEvent) => {
            if (!scatterKeyboardActiveRef.current) return;
            if (event.key === "Escape") {
                if (highlightedIndices.length === 0) return;
                event.preventDefault();
                event.stopPropagation();
                highlightRows([]);
                return;
            }
            if (event.key !== " " && event.key !== "Spacebar") return;
            event.preventDefault();
            event.stopPropagation();
            spaceHighlightModeRef.current = event.shiftKey ? "remove" : "add";
            if (spaceHighlightActiveRef.current) return;
            spaceHighlightActiveRef.current = true;
            lastSpaceHighlightedRowRef.current = -1;
            paintHoveredRow();
        },
        [highlightedIndices.length, highlightRows, paintHoveredRow],
    );
    const onScatterGlobalKeyUp = useCallback((event: KeyboardEvent) => {
        if (event.key !== " " && event.key !== "Spacebar") return;
        spaceHighlightActiveRef.current = false;
        lastSpaceHighlightedRowRef.current = -1;
    }, []);
    const clearScatterKeyboardState = useCallback(() => {
        spaceHighlightActiveRef.current = false;
        spaceHighlightModeRef.current = "add";
        lastSpaceHighlightedRowRef.current = -1;
    }, []);
    useEffect(() => {
        window.addEventListener("keydown", onScatterGlobalKeyDown);
        window.addEventListener("keyup", onScatterGlobalKeyUp);
        window.addEventListener("blur", clearScatterKeyboardState);
        return () => {
            window.removeEventListener("keydown", onScatterGlobalKeyDown);
            window.removeEventListener("keyup", onScatterGlobalKeyUp);
            window.removeEventListener("blur", clearScatterKeyboardState);
        };
    }, [onScatterGlobalKeyDown, onScatterGlobalKeyUp, clearScatterKeyboardState]);
    const setScatterKeyboardActive = useCallback(
        (active: boolean) => {
            scatterKeyboardActiveRef.current = active;
            if (!active) clearScatterKeyboardState();
        },
        [clearScatterKeyboardState],
    );
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
                // should this also refer to something other than scatterplotLayer?
                // in DeckScatterComponent we now use deckRef.current.deck.getViewports()[0]
                const p = scatterplotLayer.unproject([x, y]) as P;
                //still need to reason better about transforms...
                // const scale = 1; //this was only right when Pixels.PhysicalSizeX === regions.scale
                // const p2 = p.map(v => v * scale) as P;
                const m = modelMatrix.invert();
                const p3 = m.transform(p) as P;
                m.invert();
                return p3;
            } catch {
                //console.error("unproject error", e);
                return [0, 0];
            }
        },
        [scatterplotLayer, modelMatrix, boundingClientRect],
    );
    return useMemo(
        () => ({
            scatterplotLayer,
            greyScatterplotLayer,
            unselectedScatterplotLayer,
            densityLayer,
            getTooltip,
            modelMatrix,
            viewState,
            unproject,
            setScatterKeyboardActive,
        }),
        [
            scatterplotLayer,
            greyScatterplotLayer,
            unselectedScatterplotLayer,
            densityLayer,
            getTooltip,
            modelMatrix,
            viewState,
            unproject,
            setScatterKeyboardActive,
        ],
    );
}
