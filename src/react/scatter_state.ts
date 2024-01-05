import { PickingInfo, ScatterplotLayer } from "deck.gl/typed";
import { ScatterPlotConfig, VivRoiConfig } from "./components/VivMDVReact";
import { useChart, useDataStore, useOmeTiff } from "./context";
import { useChartID, useConfig, useParamColumns } from "./hooks";
import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState } from "react";
import { getVivId } from "./components/avivatorish/MDVivViewer";

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
    useLayoutEffect(() => {
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


/** this should be able to deal with understanding in a given context, how to scale each axis...
 * right now, hacking together something that works for adenoma/carcinoma data.
 * We should be more careful about checking that the various scale units are consistent etc.
 * In future we may want to return a matrix, and allow the user to manipulate it.
 */
export function useRegionScale() {
    const ome = useOmeTiff(); //this hook should be valid when we don't have an ome tiff as well
    const chart = useChart();
    const regionScale = chart.dataStore.regions.scale;
    //see also getPhysicalScalingMatrix
    //- consider state, matrices for image, scatterplot/other layers, and options to manipulate them
    //MDVProject.set_region_scale assumes that all regions have the same scale.
    //in the current data-set, that is true, but it could be false in the future.
    //A-priori it would make sensee for x/y columns to have metadata about their scale, but that's not the case.
    //Indeed, any type of numerical column should probably have metadata about its scale, 
    //as well as information like comments what it represents.
    //and given the other issues it could be inconsistent anyway?
    const scale = ome.metadata.Pixels.PhysicalSizeX / regionScale;
    return scale;
}
type Tooltip = (PickingInfo) => string;
export function useScatterplotLayer(): [ScatterplotLayer, Tooltip] {
    const id = useChartID();
    const chart = useChart();
    const colorBy = (chart as any).colorBy;
    const config = useConfig<ScatterPlotConfig>();
    const scale = useRegionScale();

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
        const hoverIndex = i;
        const index = data[hoverIndex];
        const valueIndex = tooltipCol.data[index];
        const value = tooltipCol.values[valueIndex];
        // return JSON.stringify({
        //     hoverIndex,
        //     index,
        //     valueIndex,
        //     value,
        // }, null, 2);
        return value;
    }, [tooltipCol, tooltipCol.data, tooltipCol.values, data]);
    const getTooltip = useCallback(
        //todo nicer tooltip interface (and review how this hook works)
        () => {
            const hoverInfo = hoverInfoRef.current;
            return hoverInfo && hoverInfo.index !== -1 && tooltipCol && `${getTooltipVal(hoverInfo.index)}`;
        },
    [hoverInfoRef, tooltipCol]);

    // debugging...
    // useEffect(() => {
    //     // get a unique set of values referred to by the tooltip column as per filtered data
    //     if (!tooltipCol) return;
    //     const indices = data.map(fi => tooltipCol.data[fi]);
    //     const indexCounts = indices.reduce((acc, i) => { acc[i] = (acc[i] ?? 0) + 1; return acc; }, {});
    //     const strings = Array.from(new Set(indices)).map(i => tooltipCol.values[i]);
    //     const stringCounts = tooltipCol.values.map(s => `${s}: ${
    //         Math.round(100*indexCounts[tooltipCol.values.indexOf(s)]/data.length)
    //     }%`);
    //     console.table(stringCounts);
    // }, [tooltipCol, tooltipCol.data, tooltipCol.values, data])

    const scatterplotLayer = useMemo(() => new ScatterplotLayer({
        id: `scatter_${getVivId(id + 'detail-react')}`, // should satisfy VivViewer, could make this tidier
        data,
        opacity,
        radiusScale: radius / scale,
        getFillColor: colorBy ?? [0, 200, 200],
        getRadius: 1,
        getPosition: (i, { target }) => {
            target[0] = cx.data[i] / scale;
            target[1] = cy.data[i] / scale;
            target[2] = 0;
            return target as unknown as Float32Array; // deck.gl types are wrong AFAICT
        },
        updateTriggers: {
            getFillColor: colorBy,
        },
        pickable: true,
        onHover: (info) => {
            hoverInfoRef.current = info;
        }        
    }), [id, data, opacity, radius, colorBy, cx, cy]);
    return [scatterplotLayer, getTooltip];
}