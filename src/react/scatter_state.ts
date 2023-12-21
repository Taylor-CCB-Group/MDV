import { PickingInfo, ScatterplotLayer } from "deck.gl/typed";
import { ScatterPlotConfig } from "./components/VivMDVReact";
import { useChart, useOmeTiff } from "./context";
import { useChartID, useConfig, useFilteredIndices, useParamColumns } from "./hooks";
import { useMemo, useState } from "react";

/** this should be able to deal with understanding in a given context, how to scale each axis...
 * right now, hacking together something that works for adenoma/carcinoma data.
 * We should be more careful about checking that the various scale units are consistent etc.
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

export function useScatterplotLayer(): [ScatterplotLayer, PickingInfo] {
    const id = useChartID();
    const colorBy = (useChart() as any).colorBy;
    const config = useConfig<ScatterPlotConfig>();
    const scale = useRegionScale();

    // seem to be reacting fine to changes, why did I think I needed to use extra autorun or reaction?
    const { opacity, radius } = config;

    const data = useFilteredIndices();
    const [cx, cy] = useParamColumns();
    const [hoverInfo, setHoverInfo] = useState<PickingInfo>(null);
    const scatterplotLayer = useMemo(() => new ScatterplotLayer({
        id: id + 'scatter-react',
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
            // todo look up text as per config
            setHoverInfo(info);
        }        
    }), [id, data, opacity, radius, colorBy, cx, cy]);
    return [scatterplotLayer, hoverInfo];
}