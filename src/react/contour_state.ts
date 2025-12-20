import type { CategoricalDataType, DataColumn, LoadedDataColumn, FieldName } from "@/charts/charts";
import { useMemo } from "react";
import {
    useCategoryFilterIndices,
    useConfig,
    useFilteredIndices,
    useFieldSpec,
    useParamColumns,
    useFieldSpecs,
} from "./hooks";
import { useDataStore } from "./context";
import { useDebounce } from "use-debounce";
import { useViewState } from "./deck_state";
import { g, isArray, toArray } from "@/lib/utils";
import { observable, autorun } from "mobx";
import type { BaseConfig } from "@/charts/BaseChart";
import type BaseChart from "@/charts/BaseChart";
import type { FieldSpec, FieldSpecs } from "@/lib/columnTypeHelpers";
import { oklch2rgb } from "@/utilities/oklch2rgb";
import { getFieldColor } from "./fieldColorManager";
import type { FieldLegendItem } from "./components/FieldContourLegend";
// import { DataFilterExtension } from '@deck.gl/extensions';

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
    fillThreshold: number;
};
export type FieldContourProps = {
    id: string;
    fill: boolean;
    bandwidth: number;
    intensity: number;
    opacity: number;
    fillThreshold: number;
    fields?: LoadedDataColumn<"double">[];
    hoveredFieldId?: FieldName | null;
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
            contourParameter ? ds.getColumnColors(contourParameter.field, {
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
    const { id, parameter, category, fill, bandwidth, intensity, opacity, fillThreshold } =
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
    // what if we don't have viv stores?
    const { zoom } = useViewState();
    // we can compensate so that we don't have radiusPixels, but it makes it very slow...
    //won't be necessary when we implement heatmap differently
    const [debounceZoom] = useDebounce(zoom, 500);

    return useMemo(() => {
        if (!category || !contourParameter) return undefined;
        //If I return a layer here rather than props, will it behave as expected?
        //not really - we want to pass this into getSublayerProps() so the id is used correctly
        const radiusPixels = 30 * bandwidth * 2 ** debounceZoom;
        // console.log('radiusPixels', radiusPixels);
        // there is an issue of the scaling of these layers e.g. with images that have been resized...
        // what is different about how we scale these layers vs other scatterplot layer?
        return {
            id,
            data,
            fillOpacity: intensity,
            contourOpacity: opacity,
            // when fill is disabled, this is some arbitrary large value, otherwise use the tweakable threshold
            contourFill: fill ? fillThreshold : 10000,
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
        contourParameter,
        intensity,
        cx,
        cy,
        colorRange,
        debounceZoom,
        bandwidth,
        fill,
        opacity,
        fillThreshold,
    ]);
}
/** pending better definition */
export type ContourLayerProps = ReturnType<typeof useCategoryContour>;
export function useFieldContour(props: FieldContourProps) {
    const { id, fill, bandwidth, intensity, opacity, fillThreshold, fields, hoveredFieldId } =
        props;
    // there's a possiblity that in future different layers of the same chart might draw from different data sources...
    // so encapsulating things like getPosition might be useful.
    const [cx, cy] = useParamColumns();
    const data = useFilteredIndices();
    const { zoom } = useViewState();
    // we can compensate so that we don't have radiusPixels, but it makes it very slow...
    //won't be necessary when we implement heatmap differently
    const [debounceZoom] = useDebounce(zoom, 500);

    return useMemo(() => {
        if (!fields) return [];
        //If I return a layer here rather than props, will it behave as expected?
        //not really - we want to pass this into getSublayerProps() so the id is used correctly
        const radiusPixels = 30 * bandwidth * 2 ** debounceZoom;
        // console.log('radiusPixels', radiusPixels);
        // there is an issue of the scaling of these layers e.g. with images that have been resized...
        // what is different about how we scale these layers vs other scatterplot layer?
        //const fieldStats = fields.reduce((field) => { ... });
        return fields.map(({ name, field: fieldId, data: fieldData, minMax }) => {
            // Adjust opacity based on hover state
            const isHovered = hoveredFieldId === fieldId;
            const hasHover = hoveredFieldId !== null && hoveredFieldId !== undefined;
            
            // Hovered field: increase opacity (clamped to 1.0)
            // Non-hovered fields: reduce opacity when another field is hovered
            const adjustedOpacity = isHovered 
                ? Math.min(1.0, opacity * 1.5)
                : hasHover 
                    ? opacity * 0.8
                    : opacity;
            
            const adjustedIntensity = isHovered 
                ? Math.min(1.0, intensity * 1.5)
                : hasHover 
                    ? intensity * 0.8
                    : intensity;
            
            return {
            id: `${id}_${name}`,//if we base this id on index rather than name we might do some transitions
            data, //todo filter sparse data
            fillOpacity: adjustedIntensity,
            contourOpacity: adjustedOpacity,
            // when fill is disabled, this is some arbitrary large value, otherwise use the tweakable threshold
            contourFill: fill ? fillThreshold : 10000,
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
                //things to consider: in some circumstances normalise range across all fields
                const [min, max] = minMax;
                const range = max - min;
                if (range === 0) return 0;
                const value = fieldData[i];
                if (value < min) return 0;
                if (value > max) return 1;
                const normalizedValue = (value - min) / range;
                return normalizedValue;
            },
            // Stable color assignment based on field identifier, not index
            colorRange: [getFieldColor(fieldId)],
            //! this scale does not adapt well to the data, and it would be nice to have meaningful units...
            radiusPixels: radiusPixels,
            debounce: 1000,
            weightsTextureSize: 128, //there could be a performance related parameter to tweak
            pickable: false,
            transitions: {
                // nb - failed to get this working
                // will need changes to the HeatmapContourExtension implementation
                // (updateState vs draw, props vs state) but attempts thus far have been fruitless
                // plan to implement that differently anyway.
                // fillOpacity: 500,
                // contourOpacity: 500
                /// getWeight transition not working either
                getWeight: {
                    duration: 1000,
                    // easing: t => t,
                },
            },
            updateTriggers: {
                getWeight: [fieldData],
                // getFilterValue: [fieldData]
            },
            // not working?
            // getFilterValue: (i: number) => fieldData[i] === 0,
            // extensions: [new DataFilterExtension({filterSize: 1})]
        };
        });
    }, [
        id,
        data,
        intensity,
        cx,
        cy,
        debounceZoom,
        bandwidth,
        fill,
        opacity,
        fillThreshold,
        fields,
        hoveredFieldId,
    ]);
}

export type ContourVisualConfig = {
    contour_fill: boolean;
    contour_fillThreshold: number;
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
    densityFields?: FieldSpecs; // don't have a way of specifying datatype here
    field_legend?: {
        display?: boolean;
        pos?: [number, number];
    };
} & ContourVisualConfig;

export function getDensitySettings(c: DualContourLegacyConfig & BaseConfig, chart: BaseChart<any>) {
    const { dataStore } = chart;
    // make it so that if we change the parameter, we get the new values in the dropdowns
    // empty array will be replaced with the new values
    const catsValues = observable.array([[] as { t: string }[], "t", "t"]) as unknown as [{ t: string }[], "t", "t"];
    
    // Create autorun that will be disposed when the settings dialog closes
    // The disposer is stored in the _disposers array on the returned GuiSpec
    // so it gets cleaned up properly when the dialog closes
    const disposer = autorun(() => {
        if (typeof c.contourParameter !== "string") {
            // as of now, categorical parameter like this is expected to be a string here
            // we would like to be operating on a version of state that just had a column object
            // complete with values, regardless of provenance, and not need to refer to dataStore
            console.error("unexpected type for contourParameter");
            return;
        }
        // getColumnValues() may throw or return undefined, especially with linked data...
        // actually, there is another issue when we don't really have a categorical column...
        try {
            const ocats = c.contourParameter ? dataStore.getColumnValues(c.contourParameter).slice() : [];
            const cats = ocats.map((x) => ({ t: x }));
            catsValues[0] = cats;
        } catch (e) {
            console.error(`error updating contour values with '${c.contourParameter}' (${e})`);
        }
    });
    
    // If we have a "category_selection" widget and it properly observed `c.contourParameter`,
    // we wouldn't need this autorun here.
    // this is related to existing selection dialog widget - both should be able to understand multitext better.
    const folderSpec = g({
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
                        //^^ maybe get Furquan to work on this.
                        label: "Contour parameter",
                        current_value: c.contourParameter || "",
                        columnType: "text",
                        func: (x) => {
                            if (x === c.contourParameter) return;
                            if (!isArray(c.param)) throw "expected param array";
                            c.contourParameter = x;
                            //ru-roh, we're not calling the 'func's... mostly we just care about reacting to the change...
                            //but setting things on config doesn't work anyway, because the dialog is based on this settings object...
                            //nb - when we switch back to a contourParameter that had categories associated, the GUI state is changing
                            // back to the old setting (but not updating the state).
                            // Shouldn't this line mean that it forgets the old categories?
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
                    g({
                        type: "multicolumn",
                        label: "Density Fields",
                        //@ts-expect-error - pending optional columns
                        current_value: c.densityFields,
                        columnType: "double",
                        func: (x) => {
                            c.densityFields = x;
                        },
                    })
                ],
            }),
            ...getContourVisualSettings(c),
            g({
                type: "folder",
                label: "Legend",
                current_value: [
                    g({
                        type: "check",
                        label: "Show Field Legend",
                        current_value: c.field_legend?.display ?? true,
                        func: (x) => {
                            if (!c.field_legend) {
                                c.field_legend = {};
                            }
                            c.field_legend.display = x;
                        },
                    }),
                ],
            }),
        ],
    });
    
    // Attach the disposer to the spec so it gets cleaned up when the settings dialog closes
    folderSpec._disposers = [disposer];
    return folderSpec;
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
            max: 50,
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
            max: 5,
            min: 0.1,
            current_value: c.contour_fillThreshold,
            continuous: true,
            label: "Fill Threshold",
            func(x) {
                c.contour_fillThreshold = x;
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
export function useLegacyDualContour(hoveredFieldId?: FieldName | null): ContourLayerProps[] {
    const config = useConfig<DualContourLegacyConfig>();
    // todo: this is currently short-circuiting for non-viv deck scatter...
    // breaking rule of hooks etc, should be fixed
    const commonProps = {
        parameter: config.contourParameter,
        fill: config.contour_fill,
        bandwidth: config.contour_bandwidth || 10,
        intensity: config.contour_intensity || 0.1,
        opacity: config.contour_opacity || 0.2,
        fillThreshold: config.contour_fillThreshold || 2,
    };
    const fields = useFieldSpecs(config.densityFields);
    const fieldContours = useFieldContour({
        ...commonProps,
        id: "fieldContours",
        fields: fields.filter(field => field.datatype === "double") as LoadedDataColumn<"double">[],
        hoveredFieldId,
    });
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

/**
 * Hook that provides legend data for field contours.
 * Returns an array of field information with their stable colors for display in a legend.
 * 
 * @param densityFields - The densityFields config from DualContourLegacyConfig
 * @returns Array of FieldLegendItem objects with name, field identifier, and color
 */
export function useFieldContourLegend(densityFields?: FieldSpecs): FieldLegendItem[] {
    const fields = useFieldSpecs(densityFields);
    
    return useMemo(() => {
        // Filter to only double fields (matching the filter in useLegacyDualContour)
        const doubleFields = fields.filter(field => field.datatype === "double");
        
        return doubleFields.map(field => ({
            name: field.name,
            field: field.field,
            color: getFieldColor(field.field),
        }));
    }, [fields]);
}
