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

function useColorRange(contourParameter: DataColumn<any>, category: string) {
    const ds = useDataStore();
    const columnColors = useMemo(() => ds.getColumnColors(contourParameter.name, {asArray: true, useValue: true}), [ds, contourParameter]);
    const categoryValueIndex = useMemo(() => {
        if (!contourParameter || !contourParameter.values) return -1;
        return contourParameter.values.indexOf(category);
    }, [contourParameter, category]);
    const colorRange = useMemo(() => {
        if (categoryValueIndex === -1) return [[0, 0, 0, 0]];
        const color = columnColors[categoryValueIndex];
        // return [[...color, 255], [...color, 255]];
        console.log('color', color, category);
        return [color]
    }, [categoryValueIndex, columnColors, category]);
    return colorRange;
}

function useCategoryFilterIndices(contourParameter: DataColumn<any>, category: string) {
    const data = useFilteredIndices();
    const categoryValueIndex = useMemo(() => {
        if (!contourParameter || !contourParameter.values) return -1;
        return contourParameter.values.indexOf(category);
    }, [contourParameter, category]);
    const filteredIndices = useMemo(() => {
        if (categoryValueIndex === -1) return [];
        return data.filter(i => contourParameter.data[i] === categoryValueIndex);
    }, [data, categoryValueIndex, contourParameter]);
    return filteredIndices;
}

export function useContour(props: ContourProps) {
    const {id, parameter, category, fill, bandwidth, intensity, opacity} = props;
    // there's a possiblity that in future different layers of the same chart might draw from different data sources...
    // so encapsulating things like getPosition might be useful.
    const [cx, cy, contourParameter] = useParamColumns();
    const data = useCategoryFilterIndices(contourParameter, category);
    // const getWeight = useContourWeight(contourParameter, category);
    const colorRange = useColorRange(contourParameter, category);

    return useMemo(() => {
        if (!category) return undefined;
        //todo if I return a layer here rather than props, will it behave as expected?
        return {
            id,
            data,
            opacity: intensity,
            getPosition: (i: number, { target }: { target: number[] | Float32Array }) => {
                target[0] = cx.data[i];
                target[1] = cy.data[i];
                target[2] = 0;
                return target;
            },
            colorRange,
        };
    }, [id, data, category, intensity, cx, cy, colorRange]);
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