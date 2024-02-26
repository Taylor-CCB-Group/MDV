import { Matrix4 } from '@math.gl/core';
import { PickingInfo } from "deck.gl/typed";
import { ScatterPlotConfig, VivRoiConfig } from "./components/VivMDVReact";
import { useChart, useDataStore } from "./context";
import { useChartID, useChartSize, useConfig, useParamColumns } from "./hooks";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { getVivId } from "./components/avivatorish/MDVivViewer";
import { OME_TIFF, useLoader, useMetadata } from "./components/avivatorish/state";
import { ViewState } from './components/VivScatterComponent';
import { getImageSize } from '@hms-dbmi/viv';
import { ScatterplotExLayer } from '../webgl/ImageArrayDeckExtension';
import { ScatterDeckExtension } from '../webgl/ScatterDeckExtension';

/**
 * Get a {Uint32Array} of the currently filtered indices.
 * When the selection changes, this will asynchronously update.
 * All users of the same data store share a reference to the same array.
 */
export function useFilteredIndices() {
    // in the case of region data, it should be filtered by that as well...
    // I really want to sort out how I use types here...
    const config = useConfig<VivRoiConfig>();
    const filterColumn = config.background_filter?.column;
    const dataStore = useDataStore();
    const [filteredIndices, setFilteredIndices] = useState(new Uint32Array());
    // todo useQuery could be better?
    useEffect(() => {
        // return
        let cancelled = false;
        if (!filterColumn) return;
        const indexPromise = dataStore.getFilteredIndices();
        const colPromise = window.mdv.chartManager?._getColumnsAsync(dataStore.name, [filterColumn]);
        Promise.all([indexPromise, colPromise]).then(([indices]) => {
        // indexPromise.then((indices) => {
            if (cancelled) return;
            if (filterColumn) {
                const col = dataStore.columnIndex[filterColumn];
                const filterValue = config.background_filter?.category;
                if (filterValue) {
                    const filterIndex = col.values.indexOf(filterValue);
                    try {
                        const filteredIndices = indices.filter(i => col.data[i] === filterIndex);
                        setFilteredIndices(filteredIndices);
                    } catch (e) {
                        console.error('error filtering indices', e);
                        return;
                    }
                    return;
                }
            }
            setFilteredIndices(indices);
        });
        // should I have a cleanup function to cancel the promise if it's not resolved
        // by the time the effect is triggered again?
        return () => {
            // if (!finished) console.log('filtered indices promise cancelled');
            cancelled = true;
        }

        // using _filteredIndicesPromise as a dependency is working reasonably well,
        // but possibly needs a bit more thought.
    }, [dataStore._filteredIndicesPromise, filterColumn, config.background_filter]);
    return filteredIndices;
}

export function useRegionScale() {
    const metadata = useMetadata();
    const chart = useChart();
    const regionScale = chart.dataStore.regions.scale;
    //see also getPhysicalScalingMatrix
    //- consider state, matrices for image, scatterplot/other layers, and options to manipulate them
    //MDVProject.set_region_scale assumes that all regions have the same scale?
    if (!metadata) return 1;
    const scale = metadata.Pixels.PhysicalSizeX / regionScale;
    return scale;
}

/** for this to be more useful as a hook will depend on state/context... */
export function useScatterModelMatrix() {
    const scale = useRegionScale();
    const s = 1/scale;
    const [modelMatrix, setModelMatrix] = useState(new Matrix4().scale(s));
    // useEffect(() => {
    //     const m = new Matrix4().scale(s);
    //     setModelMatrix(m);
    // }, [scale]);
    return {modelMatrix, setModelMatrix};
}

type Tooltip = (PickingInfo) => string;
export function useScatterplotLayer() {
    const id = useChartID();
    const chart = useChart();
    const colorBy = (chart as any).colorBy;
    const config = useConfig<ScatterPlotConfig>();

    const { opacity, radius } = config;

    const data = useFilteredIndices();
    const [cx, cy] = useParamColumns();
    const hoverInfoRef = useRef<PickingInfo>(null);
    const [highlightedObjectIndex, setHighlightedObjectIndex] = useState(-1);
    const getLineWidth = useCallback((i: number) => {
        return i === data[highlightedObjectIndex] ? 0.2*radius/scale : 0;
    }, [radius, highlightedObjectIndex, data]);

    const tooltipCol = useMemo(() => {
        if (!config.tooltip) return undefined;
        return chart.dataStore.columnIndex[config.tooltip.column]
    }, [config.tooltip.column]);
    const getTooltipVal = useCallback((i: number) => {
        if (!tooltipCol) return '';
        // careful now...
        const valueIndex = tooltipCol.data[data[i]];
        if (!tooltipCol.values) return valueIndex;
        return tooltipCol.values[valueIndex];
    }, [tooltipCol, tooltipCol?.data, tooltipCol?.values, data]);
    const getTooltip = useCallback(
        //todo nicer tooltip interface (and review how this hook works)
        () => {
            if (!config.tooltip.show) return;
            const hoverInfo = hoverInfoRef.current;
            return hoverInfo && hoverInfo.index !== -1 && `${config.tooltip.column}: ${getTooltipVal(hoverInfo.index)}`;
        },
    [hoverInfoRef, getTooltipVal, config.tooltip.show]);

    const scale = useRegionScale();
    const {modelMatrix, setModelMatrix} = useScatterModelMatrix();
    const modelMatrixRef = useRef(modelMatrix);
    // we don't want to use this if we're not using an OME_TIFF - but for now, pending better reasoning about
    // what type of chart we're in, we'll use it anyway for scaling the viewState
    const loader = useLoader() as OME_TIFF['data']; // could be refactored to have a useImageSize hook maybe
    const [chartWidth, chartHeight] = useChartSize(); //not sure we want this, potentially re-rendering too often...
    // not using as dependency for scaling viewState to data - we don't want to zoom as chart size changes
    // (at least for now - may consider making this configurable / testing it out)

    const [viewState, setViewState] = useState<ViewState>(null); //we should consider how this interacts with viv ViewerStore
    useEffect(() => {
        if (!config.zoom_on_filter) return;
        if (data.length === 0) {
            setViewState(null);
            return;
        }
        const source = Array.isArray(loader) ? loader[0] : loader;
        const { width: pixelWidth, height: pixelHeight } = getImageSize(source);
        const trueImageWidth = pixelWidth * scale;
        const trueImageHeight = pixelHeight * scale;
        // Step 1: Calculate the bounding box of the data
        // with a loop to avoid spread operator & stack overflow etc
        let minX = Number.POSITIVE_INFINITY;
        let maxX = Number.NEGATIVE_INFINITY;
        let minY = Number.POSITIVE_INFINITY;
        let maxY = Number.NEGATIVE_INFINITY;
        for (let i = 0; i < data.length; i++) {
            const d = data[i];
            // todo proper use of transform matrices
            try {
                const x = cx.data[d] / scale;
                const y = cy.data[d] / scale;
                if (x === undefined || y === undefined) {
                    console.warn('undefined data in scatterplot');
                    continue;
                }
                if (x < minX) minX = x;
                if (x > maxX) maxX = x;
                if (y < minY) minY = y;
                if (y > maxY) maxY = y;
            } catch (e) {
                console.error('error calculating bounding box', e);
                return;
            }
        }
        
        // Step 2: Calculate the center of the bounding box
        // - take into account transform matrices
        // no internalState in scatterplotLayer, so we can't use project/projectPosition
        // need to figure this out properly...
        // const [x, y] = scatterplotLayer.project([x, y, 0]);
        // const p = scatterplotLayer.projectPosition([x, y, 0]);
        // const [x, y] = p;
        // quicker to do once than for every point
        [minX, minY] = modelMatrix.transformAsPoint([minX, minY, 0]);
        [maxX, maxY] = modelMatrix.transformAsPoint([maxX, maxY, 0]);
        
        const centerX = (minX + maxX) / 2;
        const centerY = (minY + maxY) / 2;

        // modelMatrix.transformPoint([centerX, centerY, 0]);

        // Step 3: Calculate the zoom level (this is a rough approximation and may need to be adjusted)
        // should consider aspect ratio of the canvas
        // should allow some padding around the edges, based on radius, and/or some config value
        const { max, min, log2 } = Math;
        const maxZoom = 5;
        const xZoom = 1/log2((maxX - minX) / trueImageWidth);
        const yZoom = 1/log2((maxY - minY) / trueImageHeight);
        console.log('zoom', xZoom, yZoom);
        const zoomBackOff = 0;
        const zoom = min(maxZoom, max(xZoom, yZoom) - zoomBackOff);
        // Step 4: Set the view state that will be picked up by the user of this hook
        // may want to use viv viewerStore - but we don't want scatterplot to depend on that
        setViewState({
            target: [centerX, centerY, 0],
            zoom: zoom,
            minZoom: -10,
            maxZoom,
            transitionDuration: 400, // Smooth transition
            transitionEasing: x => -(Math.cos(Math.PI * x) - 1) / 2, //https://easings.net/#easeInOutSine
            // transitionInterpolator: new FlyToInterpolator({speed: 1}), //applicable for MapState - latitude is required
        });
    }, [data, cx, cy, scale]);
    const { point_shape } = config;

    const extensions = useMemo(() => {
        if (point_shape === 'circle') return [];
        return [new ScatterDeckExtension()]
    }, []);
    const scatterplotLayer = useMemo(() => new ScatterplotExLayer({
        // loaders //<< this will be interesting to learn about
        id: `scatter_${getVivId(id + 'detail-react')}`, // should satisfy VivViewer, could make this tidier
        data,
        opacity,
        radiusScale: radius,
        getFillColor: colorBy ?? [0, 200, 200],
        getRadius: 1/scale,
        getPosition: (i, { target }) => {
            target[0] = cx.data[i];
            target[1] = cy.data[i];
            target[2] = 0;
            return target as unknown as Float32Array; // deck.gl types are wrong AFAICT
        },
        modelMatrix,
        updateTriggers: {
            getFillColor: colorBy, //this is working; removing it breaks the color change...
            // modelMatrix: modelMatrix, // this is not necessary, manipulating the matrix works anyway
            // getLineWith: clickIndex, // this does not work, seems to need something like a function
            getLineWidth
        },
        pickable: true,
        onHover: (info) => {
            hoverInfoRef.current = info;
        },
        stroked: data.length < 1000, //todo make this configurable, and fix issue...
        // todo figure out why lineWidth 0 still shows up, particularly when zoomed out
        // can we make it have zero opacity? Seems like lineColor is rgb, not rgba...
        // may need a layer extension to do this properly; may want that anyway for other reasons
        getLineWidth,
        //trying to set line color to same as fill, but it makes things very muddy when zoomed out
        //getLineColor: i => i === clickIndexRef.current ? [255, 255, 255] : colorBy ?? [200, 200, 200],
        getLineColor: [255, 255, 255],
        // highlightedObjectIndex, // has some undesirable effects, but could be useful when better controlled
        onClick: ({index}) => {
            setHighlightedObjectIndex(index);
            //todo properly synchronise state with data store, allow deselection
            chart.dataStore.dataHighlighted([data[index]]);
        },
        transitions: {
            // this leads to weird behaviour when filter changes, looks ok when changing colorBy
            // getFillColor: {
            //     duration: 300,
            // },
        },
        extensions
    }), [id, data, opacity, radius, colorBy, cx, cy, highlightedObjectIndex, scale, modelMatrix]);
    return {scatterplotLayer, getTooltip, modelMatrix, setModelMatrix, modelMatrixRef, viewState};
}
