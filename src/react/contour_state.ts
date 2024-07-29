import type { DataColumn } from "@/charts/charts";
import { useMemo, useCallback } from "react";
import { useConfig, useParamColumns } from "./hooks";
import { useFilteredIndices } from "./scatter_state";
import { useChart, useDataStore } from "./context";


/** need to be clearer on which prop types are for which parts of layer spec...
 * 
 */
type ContourProps = {
    /** to be used as deck.gl sublayer id */
    id: string,
    /** the column */
    parameter: string,
    category: string,
    fill: boolean,
    bandwidth: number,
    intensity: number,
    opacity: number,
}

function useContourWeight(contourParameter: DataColumn<any>, category: string) {
    //todo consider what happens when contourParameter is not categorical/text-like.
    //down the line we may consider modulating the weight etc.
    const categoryValueIndex = useMemo(() => {
        return contourParameter.values.indexOf(category);
    }, [contourParameter, category]);
    const getContourWeight = useCallback((i: number) => {
        return contourParameter.data[i] === categoryValueIndex ? 1 : 0;
    }, [contourParameter, categoryValueIndex]);
    return getContourWeight;
}

function useColorRange(contourParameter: DataColumn<any>, category: string) {
    const ds = useDataStore();
    const columnColors = useMemo(() => ds.getColumnColors(contourParameter.name, {asArray: true, useValue: true}), [ds, contourParameter]);
    const categoryValueIndex = useMemo(() => {
        return contourParameter.values.indexOf(category);
    }, [contourParameter, category]);
    const colorRange = useMemo(() => {
        const color = columnColors[categoryValueIndex];
        // return [[...color, 255], [...color, 255]];
        console.log('color', color, category);
        return [color]
    }, [categoryValueIndex, columnColors, category]);
    return colorRange;
}

export function useContour(props: ContourProps) {
    const {id, parameter, category, fill, bandwidth, intensity, opacity} = props;
    // there's a possiblity that in future different layers of the same chart might draw from different data sources...
    // so encapsulating things like getPosition might be useful.
    const data = useFilteredIndices();
    const [cx, cy, contourParameter] = useParamColumns();
    const getWeight = useContourWeight(contourParameter, category);
    const colorRange = useColorRange(contourParameter, category);

    return useMemo(() => ({
        id,
        data,
        getWeight,
        opacity: intensity,
        getPosition: (i: number, { target }: { target: number[] | Float32Array }) => {
            target[0] = cx.data[i];
            target[1] = cy.data[i];
            target[2] = 0;
            return target;
        },
        colorRange,
        updateTriggers: {
            getWeight
        },
        // debounceTimeout: 1000,
    }), [id, data, getWeight, intensity, cx, cy, colorRange]);
}

/** In future I think we want something more flexible & expressive,
 * but this should be somewhat compatible with the previous implementation 
 */
export type DualContourLegacyConfig = {
    contourParameter?: string, //this is param[2] in the original code
    category1?: string,
    category2?: string,
    contour_fill: boolean,
    contour_bandwidth: number,
    contour_intensity: number,
    contour_opacity: number,
}


/** */
export function useLegacyDualContour(){
    const config = useConfig<DualContourLegacyConfig>();
    const commonProps = {
        parameter: config.contourParameter,
        fill: config.contour_fill,
        bandwidth: config.contour_bandwidth,
        intensity: config.contour_intensity,
        opacity: config.contour_opacity,
    };
    const contour1 = useContour({ ...commonProps, id: 'contour1', category: config.category1 });
    const contour2 = useContour({ ...commonProps, id: 'contour2', category: config.category2 });
    const stableArray = useMemo(() => [contour1, contour2], [contour1, contour2]);
    return stableArray;
}