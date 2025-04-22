import type { CategoricalDataType, DataColumn, LoadedDataColumn } from "@/charts/charts";
import { useMemo } from "react";
import {
    useCategoryFilterIndices,
    useConfig,
    useFilteredIndices,
    useFieldSpec,
    useParamColumns,
} from "./hooks";
import { useDataStore } from "./context";
import { useDebounce } from "use-debounce";
import { useViewState } from "./deck_state";
import { g, isArray, toArray } from "@/lib/utils";
import { observable } from "mobx";
import type { BaseConfig } from "@/charts/BaseChart";
import type BaseChart from "@/charts/BaseChart";
import type { FieldSpec } from "@/lib/columnTypeHelpers";

/** need to be clearer on which prop types are for which parts of layer spec...
 *
 */
export type CategoryContourProps = {
    /** to be used as deck.gl sublayer id */
    id: string;
    /** the column - if not present, should we not filter the data? at the moment, we reduce to nothing
     * maybe it should not be optional, in which case we need to re-arrange how hooks are called
     */
    parameter?: FieldSpec;
    category?: string | string[];
    fill: boolean;
    bandwidth: number;
    intensity: number;
    opacity: number;
};
export type FieldContourProps = {
    id: string;
    fill: boolean;
    bandwidth: number;
    intensity: number;
    opacity: number;
    fields: LoadedDataColumn<"double">[];
}
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

export function useCategoryContour(props: CategoryContourProps) {
    const { id, parameter, category, fill, bandwidth, intensity, opacity } =
        props;
    // there's a possiblity that in future different layers of the same chart might draw from different data sources...
    // so encapsulating things like getPosition might be useful.
    const [cx, cy] = useParamColumns();
    // this has a regression... while in the process of loading data, it will return undefined which we don't handle well.
    const contourParameter = useFieldSpec(parameter);
    const data = useCategoryFilterIndices(contourParameter, category);
    // const getWeight = useContourWeight(contourParameter, category);
    const colorRange = useColorRange(contourParameter, category);
    // what if we don't have viv stores?
    const { zoom } = useViewState();
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
/** pending better definition */
export type ContourLayerProps = ReturnType<typeof useCategoryContour>;
export function useFieldContour(props: FieldContourProps) {
    const { id, fill, bandwidth, intensity, opacity, fields } =
        props;
    // there's a possiblity that in future different layers of the same chart might draw from different data sources...
    // so encapsulating things like getPosition might be useful.
    const [cx, cy] = useParamColumns();
    const data = useFilteredIndices();
    // const getWeight = useContourWeight(contourParameter, category);
    const colorRange = viridis;
    const { zoom } = useViewState();
    // we can compensate so that we don't have radiusPixels, but it makes it very slow...
    //won't be necessary when we implement heatmap differently
    const [debounceZoom] = useDebounce(zoom, 500);

    return useMemo(() => {
        //If I return a layer here rather than props, will it behave as expected?
        //not really - we want to pass this into getSublayerProps() so the id is used correctly
        const radiusPixels = 30 * bandwidth * 2 ** debounceZoom;
        // console.log('radiusPixels', radiusPixels);
        // there is an issue of the scaling of these layers e.g. with images that have been resized...
        // what is different about how we scale these layers vs other scatterplot layer?
        return fields.map(({ name, data: fieldData, minMax }) => ({
            id: `${id}_${name}`,//if we base this id on index rather than name we might do some transitions
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
            getWeight: (i: number) => {
                //potential for this to be animated in fun ways...
                //phase shift for each field in shader...

                //normalization pending refactor/design
                const [min, max] = minMax;
                const range = max - min;
                const value = fieldData[i];
                if (value < min) return 0;
                if (value > max) return 1;
                const normalizedValue = (value - min) / range;
                return normalizedValue;
            },
            colorRange,
            radiusPixels,
            debounce: 1000,
            weightsTextureSize: 512, //there could be a performance related parameter to tweak
            pickable: false,
            updateTriggers: {
                getWeight: [fieldData],
            }
        }));
    }, [
        id,
        data,
        intensity,
        cx,
        cy,
        colorRange,
        debounceZoom,
        bandwidth,
        fill,
        opacity,
        fields,
    ]);
}

export type ContourVisualConfig = {
    contour_fill: boolean;
    /** KDE bandwidth/radius. todo: units */
    contour_bandwidth: number;
    contour_intensity: number;
    contour_opacity: number;
}
/** In future I think we want something more flexible & expressive,
 * but this should be somewhat compatible with the previous implementation
 */
export type DualContourLegacyConfig = {
    //config vs state again... could we make this a DataColumn?
    contourParameter?: FieldSpec; //this is param[2] in the original code
    category1?: string | string[];
    category2?: string | string[];
} & ContourVisualConfig;

export function getDensitySettings(c: DualContourLegacyConfig & BaseConfig, chart: BaseChart<any>) {
    const { dataStore } = chart;
    // make it so that if we change the parameter, we get the new values in the dropdowns
    // empty array will be replaced with the new values
    const catsValues = observable.array([[] as { t: string }[], "t", "t"]) as unknown as [{ t: string }[], "t", "t"];
    // this autorun will be disposed when the chart is disposed... really it should be tied to the settings dialog
    chart.mobxAutorun(() => {
        if (typeof c.contourParameter !== "string") {
            // as of now, categorical parameter like this is expected to be a string here
            // we would like to be operating on a version of state that just had a column object
            // complete with values, regardless of provenance, and not need to refer to dataStore
            console.error("unexpected type for contourParameter");
            return;
        }
        const ocats = c.contourParameter ? dataStore.getColumnValues(c.contourParameter).slice() : [];
        const cats = ocats.map((x) => ({ t: x }));
        catsValues[0] = cats;
    });
    return g({
        type: "folder",
        label: "Density Visualisation",
        current_value: [
            g({
                type: "folder",
                label: "Category selection",
                current_value: [
                    //maybe 2-spaces format is better...
                    g({
                        type: "column", //todo, make this "column" and fix odd behaviour with showing the value...
                        //todo: make the others be "category_selection" or something (which we don't have yet as a GuiSpec type)
                        label: "Contour parameter",
                        current_value: c.contourParameter || c.param[2],
                        columnType: "text",
                        func: (x) => {
                            if (x === c.contourParameter) return;
                            if (!isArray(c.param)) throw "expected param array";
                            c.contourParameter = c.param[2] = x;
                            //ru-roh, we're not calling the 'func's... mostly we just care about reacting to the change...
                            //but setting things on config doesn't work anyway, because the dialog is based on this settings object...
                            c.category1 = c.category2 = [];
                        },
                    }),
                    g({
                        type: "multidropdown",
                        label: "Contour Category 1",
                        current_value: toArray(c.category1 || "None"),
                        values: catsValues,
                        func(x) {
                            // if (x === "None") x = null;
                            c.category1 = x;
                        },
                    }),
                    g({
                        type: "multidropdown",
                        label: "Contour Category 2",
                        current_value: toArray(c.category2 || "None"),
                        values: catsValues,
                        func(x) {
                            // if (x === "None") x = null;
                            c.category2 = x;
                        },
                    }),
                ],
            }),
            ...getContourVisualSettings(c)
        ],
    });
}

/**
 * Gets settings for tweaking visual properties of contours.
 * @param c - contour visual config, this might be a `chart.config` or in future a `layer.config`;
 * shouldn't require any changes to this function if we have a more flexible nested config
 * @returns - array of visual settings
 */
export function getContourVisualSettings(c: ContourVisualConfig) {
    return [
        g({
            type: "slider",
            max: 25,
            min: 1,

            // doc: this.__doc__, //why?
            current_value: c.contour_bandwidth,
            label: "KDE Bandwidth",
            continuous: true,
            func(x) {
                c.contour_bandwidth = x;
            },
        }),
        g({
            label: "Fill Contours",
            type: "check",
            current_value: c.contour_fill,
            func(x) {
                c.contour_fill = x;
            },
        }),
        g({
            type: "slider",
            max: 1,
            min: 0,
            current_value: c.contour_intensity,
            continuous: true,
            label: "Fill Intensity",
            func(x) {
                c.contour_intensity = x;
            },
        }),
        g({
            type: "slider",
            max: 1,
            min: 0,
            current_value: c.contour_opacity,
            continuous: false, //why so slow?
            label: "Contour opacity",
            func(x) {
                c.contour_opacity = x ** 3;
            },
        }),
    ]
}

/**
 * In future we will want more flexible array of contours.
 * Dual-contour is a special case of this - may be useful in terms
 * of how it relates to cell-pair interation links
 */
export function useLegacyDualContour(): ContourLayerProps[] {
    const config = useConfig<DualContourLegacyConfig>();
    // todo: this is currently short-circuiting for non-viv deck scatter...
    // breaking rule of hooks etc, should be fixed
    const commonProps = {
        parameter: config.contourParameter,
        fill: config.contour_fill || true,
        bandwidth: config.contour_bandwidth || 10,
        intensity: config.contour_intensity || 0.1,
        opacity: config.contour_opacity || 0.2,
    };
    const [cx, cy, ...fields] = useParamColumns() as LoadedDataColumn<"double">[];
    const fieldContours = useFieldContour({
        ...commonProps,
        id: "fieldContours",
        // the spread above is wrong particularly when param[2] is categorical
        fields: fields.filter(field => field.datatype === "double")
    });
    //@ts-expect-error Type 'SharedArrayBuffer' is missing the following properties from type 'ArrayBuffer'?
    if (!config.contourParameter) return fieldContours;
    const contour1 = useCategoryContour({
        ...commonProps,
        id: "contour1",
        category: config.category1,
    });
    const contour2 = useCategoryContour({
        ...commonProps,
        id: "contour2",
        category: config.category2,
    });
    const stableArray = useMemo(
        () => [contour1, contour2, ...fieldContours].filter(v => v !== undefined),
        [contour1, contour2, fieldContours],
    );
    //@ts-expect-error Type 'SharedArrayBuffer' is missing the following properties from type 'ArrayBuffer'?
    return stableArray;
}
