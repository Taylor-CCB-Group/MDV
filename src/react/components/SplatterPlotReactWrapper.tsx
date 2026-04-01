import DeckGL from "@deck.gl/react";
import { OrthographicView } from "@deck.gl/core";
import { HeatmapLayer } from "deck.gl";
import { observer } from "mobx-react-lite";
import BaseChart from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import type { CategoricalDataType, LoadedDataColumn, NumberDataType } from "../../charts/charts";
import type DataStore from "@/datastore/DataStore";
import { g } from "@/lib/utils";
import { useChartSize, useConfig, useFieldSpec, useFieldSpecs, useFilteredIndices, useParamColumns } from "../hooks";
import { getDensityVisualisationFolder } from "../contour_state";
import { useCallback, useMemo } from "react";
import { TriangleLayerContours } from "@/webgl/HeatmapContourExtension";
import { getFieldColor } from "../fieldColorManager";
import {
    adaptSplatterConfig,
    getSharedSplatterViewState,
    getSplatterLayout,
    getVisibleSplatterCategories,
    matchesSplatterView,
    type SplatterPlotConfig,
} from "./splatterPlotUtils";

function buildFieldColorRange(fieldId: string) {
    const color = getFieldColor(fieldId);
    return Array.from({ length: 6 }, () => [...color]);
}

const SplatterPlot = observer(function SplatterPlot() {
    const [width, height] = useChartSize();
    const [cx, cy] = useParamColumns() as LoadedDataColumn<NumberDataType>[];
    const config = useConfig<SplatterPlotConfig>();
    const filteredRows = useFilteredIndices();
    const categoryColumn = useFieldSpec(config.category) as LoadedDataColumn<CategoricalDataType> | undefined;
    const densityFields = useFieldSpecs(config.densityFields).filter(
        (field): field is LoadedDataColumn<NumberDataType> =>
            field.datatype === "double" || field.datatype === "integer" || field.datatype === "int32",
    );
    const hasConfiguredDensityFields = Array.isArray(config.densityFields) && config.densityFields.length > 0;

    const visibleCategories = useMemo(
        () => getVisibleSplatterCategories(categoryColumn, filteredRows),
        [categoryColumn, filteredRows],
    );
    const layout = useMemo(
        () =>
            getSplatterLayout(
                width,
                height,
                visibleCategories.map((category) => category.label),
                densityFields.map((field) => field.name),
            ),
        [width, height, visibleCategories, densityFields],
    );
    const sharedViewState = useMemo(
        () => getSharedSplatterViewState(cx, cy, filteredRows, layout.cellWidth, layout.cellHeight),
        [cx, cy, filteredRows, layout.cellWidth, layout.cellHeight],
    );
    const radiusPixels = useMemo(
        () => 30 * config.contour_bandwidth * 2 ** Number(sharedViewState.zoom ?? 0),
        [config.contour_bandwidth, sharedViewState.zoom],
    );

    const views = useMemo(
        () =>
            visibleCategories.flatMap((category, rowIndex) =>
                densityFields.map(
                    (field, columnIndex) =>
                        new OrthographicView({
                            id: `splatter-${category.categoryIndex}-${field.field}`,
                            width: layout.cellWidth,
                            height: layout.cellHeight,
                            x: columnIndex * layout.cellWidth,
                            y: rowIndex * layout.cellHeight,
                            flipY: false,
                        }),
                ),
            ),
        [visibleCategories, densityFields, layout.cellWidth, layout.cellHeight],
    );

    const viewState = useMemo(
        () =>
            Object.fromEntries(
                visibleCategories.flatMap((category) =>
                    densityFields.map((field) => [
                        `splatter-${category.categoryIndex}-${field.field}`,
                        { ...sharedViewState },
                    ]),
                ),
            ),
        [visibleCategories, densityFields, sharedViewState],
    );

    const layers = useMemo(
        () =>
            visibleCategories.flatMap((category) =>
                densityFields.map((field) => {
                    const [min, max] = field.minMax;
                    const range = max - min;
                    const layerId = `splatter-${category.categoryIndex}-${field.field}`;

                    return new HeatmapLayer({
                        id: layerId,
                        viewId: layerId,
                        data: category.rows,
                        fillOpacity: config.contour_intensity,
                        contourOpacity: config.contour_opacity,
                        contourFill: config.contour_fill ? config.contour_fillThreshold : 10000,
                        colorRange: buildFieldColorRange(field.field),
                        radiusPixels,
                        debounce: 250,
                        weightsTextureSize: 128,
                        pickable: false,
                        getPosition: (
                            rowIndex: number,
                            { target }: { target: number[] | Float32Array },
                        ) => {
                            target[0] = cx.data[rowIndex];
                            target[1] = cy.data[rowIndex];
                            target[2] = 0;
                            return target;
                        },
                        getWeight: (rowIndex: number) => {
                            if (range === 0) return 0;
                            const value = field.data[rowIndex];
                            if (!Number.isFinite(value) || value <= min) return 0;
                            if (value >= max) return 1;
                            return (value - min) / range;
                        },
                        updateTriggers: {
                            getPosition: [cx.data, cy.data],
                            getWeight: [field.data],
                        },
                        _subLayerProps: {
                            triangle: {
                                type: TriangleLayerContours,
                            },
                            "triangle-layer": {
                                contourOpacity: config.contour_opacity,
                                contourFill: config.contour_fill ? config.contour_fillThreshold : 10000,
                                fillOpacity: config.contour_intensity,
                            },
                        },
                    } as any);
                }),
            ),
        [
            visibleCategories,
            densityFields,
            config.contour_fill,
            config.contour_fillThreshold,
            config.contour_intensity,
            config.contour_opacity,
            radiusPixels,
            cx.data,
            cy.data,
        ],
    );
    const layerFilter = useCallback(
        ({ layer, viewport }: { layer: { id: string }; viewport: { id: string } }) =>
            matchesSplatterView(layer.id, viewport.id),
        [],
    );

    if (!config.category) {
        return <div className="flex h-full items-center justify-center text-sm">Choose a category column to build the splatter plot.</div>;
    }
    if (!categoryColumn) {
        return <div className="flex h-full items-center justify-center text-sm">Loading category data...</div>;
    }
    if (densityFields.length === 0) {
        return (
            <div className="flex h-full items-center justify-center text-sm">
                {hasConfiguredDensityFields ? "Loading density fields..." : "Choose one or more density fields to build the splatter plot."}
            </div>
        );
    }
    if (visibleCategories.length === 0) {
        return <div className="flex h-full items-center justify-center text-sm">No categories remain after the current filters.</div>;
    }

    return (
        <div className="relative h-full w-full">
            <div
                className="absolute left-0 top-0 grid"
                style={{
                    top: `${layout.headerHeight}px`,
                    width: `${layout.labelWidth}px`,
                    height: `${layout.plotHeight}px`,
                    gridTemplateRows: `repeat(${visibleCategories.length}, minmax(0, 1fr))`,
                }}
            >
                {visibleCategories.map((category) => (
                    <div
                        key={`row-label-${category.categoryIndex}`}
                        className="flex items-center justify-end pr-3 text-right text-xs"
                        title={category.label}
                    >
                        <span className="max-w-full truncate">{category.label === "" ? "none" : category.label}</span>
                    </div>
                ))}
            </div>
            <div
                className="absolute top-0 grid"
                style={{
                    left: `${layout.labelWidth}px`,
                    width: `${layout.plotWidth}px`,
                    height: `${layout.headerHeight}px`,
                    gridTemplateColumns: `repeat(${densityFields.length}, minmax(0, 1fr))`,
                }}
            >
                {densityFields.map((field) => {
                    const color = getFieldColor(field.field);
                    return (
                        <div key={`column-label-${field.field}`} className="relative overflow-visible">
                            <div
                                className="absolute bottom-2 left-2 whitespace-nowrap text-xs"
                                style={{
                                    color: `rgb(${color[0]}, ${color[1]}, ${color[2]})`,
                                    transform: "rotate(-35deg)",
                                    transformOrigin: "left bottom",
                                }}
                                title={field.name}
                            >
                                {field.name === "" ? "none" : field.name}
                            </div>
                        </div>
                    );
                })}
            </div>
            <div
                className="absolute"
                style={{
                    left: `${layout.labelWidth}px`,
                    top: `${layout.headerHeight}px`,
                    width: `${layout.plotWidth}px`,
                    height: `${layout.plotHeight}px`,
                }}
            >
                <DeckGL
                    controller={false}
                    layerFilter={layerFilter}
                    layers={layers}
                    views={views}
                    viewState={viewState as any}
                    getCursor={() => "default"}
                    useDevicePixels={true}
                />
                <div
                    className="absolute inset-0 grid pointer-events-none"
                    style={{
                        gridTemplateColumns: `repeat(${densityFields.length}, minmax(0, 1fr))`,
                        gridTemplateRows: `repeat(${visibleCategories.length}, minmax(0, 1fr))`,
                    }}
                >
                    {visibleCategories.flatMap((category) =>
                        densityFields.map((field) => (
                            <div
                                key={`grid-${category.categoryIndex}-${field.field}`}
                                style={{
                                    border: "1px solid rgba(148, 163, 184, 0.35)",
                                }}
                            />
                        )),
                    )}
                </div>
            </div>
        </div>
    );
});

class DeckSplatterReact extends BaseReactChart<SplatterPlotConfig> {
    constructor(dataStore: DataStore, div: HTMLDivElement, originalConfig: SplatterPlotConfig) {
        const config = adaptSplatterConfig(originalConfig);
        if (!config.title) {
            config.title =
                typeof config.category === "string"
                    ? (dataStore.getColumnName(config.category) ?? config.category)
                    : "Splatter Plot";
        }
        super(dataStore, div, config, SplatterPlot);
    }

    getSettings() {
        const c = this.config;
        return super.getSettings().concat([
            getDensityVisualisationFolder(c, {
                categorySelectionControls: [
                    g({
                        type: "column",
                        label: "Categories on y-axis",
                        current_value: c.category || "",
                        columnType: ["text", "text16"],
                        func: (x) => {
                            c.category = x;
                        },
                    }),
                    g({
                        type: "multicolumn",
                        label: "Density Fields",
                        current_value: c.densityFields || [],
                        columnType: "number",
                        func: (x) => {
                            c.densityFields = x;
                        },
                    }),
                ],
            }),
        ]);
    }
}

BaseChart.types["DeckSplatter"] = {
    name: "Splatter Plot",
    class: DeckSplatterReact,
    allow_user_add: true,
    required: (ds) =>
        ds.getColumnList("number").length >= 3 && ds.getColumnList(["text", "text16"]).length >= 1,
    params: [
        {
            type: "number",
            name: "x axis",
        },
        {
            type: "number",
            name: "y axis",
        },
        {
            type: "_multi_column:number",
            name: "density fields",
        },
    ],
    extra_controls: (ds) => [
        {
            type: "dropdown",
            name: "category",
            label: "Categories on y-axis",
            values: ds.getColumnList(["text", "text16"]).map((column) => ({
                name: column.name,
                value: column.field,
            })),
            defaultVal: ds.getColumnList(["text", "text16"])[0]?.field,
        },
    ],
    init: (config, _dataStore, extraControls) => {
        const [x, y, ...densityFields] = config.param || [];
        config.param = [x, y];
        config.category = extraControls.category;
        config.densityFields = densityFields;
    },
};

export default 42;
