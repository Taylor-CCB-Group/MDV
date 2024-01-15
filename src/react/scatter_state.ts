import { Matrix4 } from '@math.gl/core';
import { PickingInfo, ScatterplotLayer } from "deck.gl/typed";
import { ScatterPlotConfig, VivRoiConfig } from "./components/VivMDVReact";
import { useChart, useDataStore } from "./context";
import { useChartID, useConfig, useParamColumns } from "./hooks";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { getVivId } from "./components/avivatorish/MDVivViewer";
import { useMetadata } from "./components/avivatorish/state";

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
    useEffect(() => {
        // return
        let cancelled = false;
        let finished = false;
        const promise = dataStore.getFilteredIndices();
        promise.then((indices) => {
            if (cancelled) return;
            finished = true;
            if (filterColumn) {
                const col = dataStore.columnIndex[filterColumn];
                const filterValue = config.background_filter?.category;
                if (filterValue) {
                    const filterIndex = col.values.indexOf(filterValue);
                    const filteredIndices = indices.filter(i => col.data[i] === filterIndex);
                    setFilteredIndices(filteredIndices);
                    return;
                }
            }
            setFilteredIndices(indices);
        });
        // should I have a cleanup function to cancel the promise if it's not resolved
        // by the time the effect is triggered again?
        return () => {
            if (!finished) console.log('filtered indices promise cancelled');
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

    // seem to be reacting fine to changes, why did I think I needed to use extra autorun or reaction?
    const { opacity, radius } = config;

    const data = useFilteredIndices();
    const [cx, cy] = useParamColumns();
    const hoverInfoRef = useRef<PickingInfo>(null);

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
            const hoverInfo = hoverInfoRef.current;
            return hoverInfo && hoverInfo.index !== -1 && `${getTooltipVal(hoverInfo.index)}`;
        },
    [hoverInfoRef, getTooltipVal]);

    const scale = useRegionScale();
    const {modelMatrix, setModelMatrix} = useScatterModelMatrix();
    const modelMatrixRef = useRef(modelMatrix);
    const scatterplotLayer = useMemo(() => new ScatterplotLayer({
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
            getFillColor: colorBy,
            modelMatrix: modelMatrix
        },
        pickable: true,
        onHover: (info) => {
            hoverInfoRef.current = info;
        }        
    }), [id, data, opacity, radius, colorBy, cx, cy]);
    return {scatterplotLayer, getTooltip, modelMatrix, setModelMatrix, modelMatrixRef};
}