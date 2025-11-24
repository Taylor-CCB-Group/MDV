import type { CategoricalDataType, DataColumn } from "@/charts/charts";
import { useMemo } from "react";
import {
    useCategoryFilterIndices,
    useConfig,
    useFieldSpec,
    useParamColumns,
} from "./hooks";
import { useDataStore } from "./context";
import { useViewerStore } from "./components/avivatorish/state";
import { useDebounce } from "use-debounce";

/** need to be clearer on which prop types are for which parts of layer spec...
 *
 */
export type ContourProps = {
    /** to be used as deck.gl sublayer id */
    id: string;
    /** the column - if not present, should we not filter the data? at the moment, we reduce to nothing
     * maybe it should not be optional, in which case we need to re-arrange how hooks are called
     */
    parameter?: string;
    category?: string | string[];
    fill: boolean;
    bandwidth: number;
    intensity: number;
    opacity: number;
};

function rgb(
    r: number,
    g: number,
    b: number,
    a = 255,
): [number, number, number, number] {
    return [r, g, b, a];
}
// this is not the way to do it...
const contourColors = Array.from({ length: 200 }, (_, i) => {
    const v = i % 20 <= 1 ? 255 : 0;
    return rgb(v, v, v, v);
});
const viridis = [
    rgb(0, 47, 97),
    rgb(0, 95, 133),
    rgb(0, 139, 152),
    rgb(0, 181, 153),
    rgb(24, 220, 130),
    rgb(151, 245, 84),
    rgb(255, 255, 0),
] as const;

function useColorRange(
    contourParameter: DataColumn<CategoricalDataType> | undefined,
    category: string | string[] | undefined,
) {
    const ds = useDataStore();
    const columnColors = useMemo(
        () =>
            contourParameter ? ds.getColumnColors(contourParameter.name, {
                asArray: true,
                useValue: true,
            }) : viridis,
        [ds, contourParameter],
    );
    /** 
     * if the category refers to a specific value, 
     * return its index in `category.values` (and associated `columnColors`).
     * otherwise return `-1` indicating that we should use the default color range
     * (currently hardcoded to `viridis`)
     */
    const categoryValueIndex = useMemo(() => {
        // if (!category) return contourParameter.values.map((_, i) => i); //NO: -1 is a clue to use general 'viridis' color range atm.
        if (!contourParameter || !category) return -1;
        //we could do something different here... would need more clever color handling on the receiving end
        if (Array.isArray(category)) return category.length > 1 ? -1 : contourParameter.values.indexOf(category[0]);
        return contourParameter.values.indexOf(category);
    }, [contourParameter, category]);
    const colorRange = useMemo(() => {
        if (categoryValueIndex === -1) return viridis;
        const color = columnColors[categoryValueIndex];
        // return [[...color, 255], [...color, 255]];
        console.log("color", color, category);
        // return [color];
        // workaround for https://github.com/visgl/deck.gl/issues/9219
        // always use same length array so it doesn't delete the texture
        return new Array(viridis.length).fill(color);
    }, [categoryValueIndex, columnColors, category]);
    return colorRange;
}

// this should be moved elsewhere - category_state.ts? or hooks.ts - as long as HMR works

export function useContour(props: ContourProps) {
    const { id, parameter, category, fill, bandwidth, intensity, opacity } =
        props;
    // there's a possiblity that in future different layers of the same chart might draw from different data sources...
    // so encapsulating things like getPosition might be useful.
    const [cx, cy] = useParamColumns();
    // this has a regression... while in the process of loading data, it will return undefined which we don't handle well.
    // (could do with some more testing but I think this is less problematic now?)
    const contourParameter = useFieldSpec(parameter);
    const data = useCategoryFilterIndices(contourParameter, category);
    // const getWeight = useContourWeight(contourParameter, category);
    const colorRange = useColorRange(contourParameter, category);
    const { zoom } = useViewerStore((store) => store.viewState) ?? { zoom: 0 };
    // we can compensate so that we don't have radiusPixels, but it makes it very slow...
    //won't be necessary when we implement heatmap differently
    const [debounceZoom] = useDebounce(zoom, 500);

    return useMemo(() => {
        if (!category) return undefined;
        //If I return a layer here rather than props, will it behave as expected?
        //not really - we want to pass this into getSublayerProps() so the id is used correctly
        const radiusPixels = 30 * bandwidth * 2 ** debounceZoom;
        // console.log('radiusPixels', radiusPixels);
        // there is an issue of the scaling of these layers e.g. with images that have been resized...
        // what is different about how we scale these layers vs other scatterplot layer?
        return {
            id,
            data,
            opacity: fill ? intensity : 0,
            contourOpacity: opacity,
            getPosition: (
                i: number,
                { target }: { target: number[] | Float32Array },
            ) => {
                target[0] = cx.data[i];
                target[1] = cy.data[i];
                target[2] = 0;
                return target;
            },
            colorRange,
            radiusPixels,
            debounce: 1000,
            weightsTextureSize: 512, //there could be a performance related parameter to tweak
            pickable: false,
        };
    }, [
        id,
        data,
        category,
        intensity,
        cx,
        cy,
        colorRange,
        debounceZoom,
        bandwidth,
        fill,
        opacity,
    ]);
}

/** In future I think we want something more flexible & expressive,
 * but this should be somewhat compatible with the previous implementation
 */
export type DualContourLegacyConfig = {
    contourParameter?: string; //this is param[2] in the original code
    category1?: string | string[];
    category2?: string | string[];
    contour_fill: boolean;
    contour_bandwidth: number;
    contour_intensity: number;
    contour_opacity: number;
};

/**
 * In future we will want more flexible array of contours.
 * Dual-contour is a special case of this - may be useful in terms
 * of how it relates to cell-pair interation links
 */
export function useLegacyDualContour() {
    const config = useConfig<DualContourLegacyConfig>();
    // todo: this is currently short-circuiting for non-viv deck scatter...
    // breaking rule of hooks etc, should be fixed
    if (!config.contourParameter) return [];
    const commonProps = {
        parameter: config.contourParameter,
        fill: config.contour_fill,
        bandwidth: config.contour_bandwidth,
        intensity: config.contour_intensity,
        opacity: config.contour_opacity,
    };
    const contour1 = useContour({
        ...commonProps,
        id: "contour1",
        category: config.category1,
    });
    const contour2 = useContour({
        ...commonProps,
        id: "contour2",
        category: config.category2,
    });
    const stableArray = useMemo(
        () => [contour1, contour2],
        [contour1, contour2],
    );
    return stableArray;
}
