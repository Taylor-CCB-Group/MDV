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
import { useCallback, useMemo, useState } from "react";
import { TriangleLayerContours } from "@/webgl/HeatmapContourExtension";
import { getFieldColor } from "../fieldColorManager";
import { useHighlightCategoryRows } from "../selectionHooks";
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

function getCellKey(rowIndex: number, columnIndex: number) {
    return `${rowIndex}:${columnIndex}`;
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
    const [hoveredCell, setHoveredCell] = useState<{ rowIndex: number; columnIndex: number } | null>(null);
    const [highlightedCells, setHighlightedCells] = useState<Set<string>>(() => new Set());

    const visibleCategories = useMemo(
        () => getVisibleSplatterCategories(categoryColumn, filteredRows),
        [categoryColumn, filteredRows],
    );
    const highlightCategoryRows = useHighlightCategoryRows(visibleCategories);
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
    const updateHighlightedCells = useCallback((rowIndex: number, columnIndex: number, toggle = false) => {
        const cellKey = getCellKey(rowIndex, columnIndex);
        setHighlightedCells((current) => {
            const next = toggle ? new Set(current) : new Set<string>();
            if (toggle) {
                if (next.has(cellKey)) next.delete(cellKey);
                else next.add(cellKey);
            } else {
                next.add(cellKey);
            }

            const highlightedCategoryIndices = Array.from(next, (selectedCell) =>
                Number.parseInt(selectedCell.split(":")[0] ?? "-1", 10),
            ).filter((index) => index >= 0);
            highlightCategoryRows(highlightedCategoryIndices);
            return next;
        });
    }, [highlightCategoryRows]);
    const clearHighlightedCells = useCallback(() => {
        setHighlightedCells((current) => (current.size === 0 ? current : new Set()));
        highlightCategoryRows([]);
    }, [highlightCategoryRows]);
    const isCellHighlighted = useCallback(
        (rowIndex: number, columnIndex: number) => highlightedCells.has(getCellKey(rowIndex, columnIndex)),
        [highlightedCells],
    );
    const hoveredRowIndex = hoveredCell?.rowIndex ?? null;
    const hoveredColumnIndex = hoveredCell?.columnIndex ?? null;
    const highlightedRowIndices = useMemo(
        () =>
            new Set(
                Array.from(highlightedCells, (cellKey) =>
                    Number.parseInt(cellKey.split(":")[0] ?? "-1", 10),
                ).filter((index) => index >= 0),
            ),
        [highlightedCells],
    );
    const highlightedColumnIndices = useMemo(
        () =>
            new Set(
                Array.from(highlightedCells, (cellKey) =>
                    Number.parseInt(cellKey.split(":")[1] ?? "-1", 10),
                ).filter((index) => index >= 0),
            ),
        [highlightedCells],
    );
    const hoveredCategory = hoveredRowIndex === null ? null : visibleCategories[hoveredRowIndex];
    const hoveredField = hoveredColumnIndex === null ? null : densityFields[hoveredColumnIndex];
    const gridLineColor = "rgba(148, 163, 184, 0.16)";

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
        <div
            className="relative h-full w-full"
            onKeyDown={(event) => {
                if (event.key !== "Escape") return;
                event.preventDefault();
                event.stopPropagation();
                clearHighlightedCells();
            }}
        >
            <div
                className="absolute left-0 top-0 grid"
                style={{
                    top: `${layout.headerHeight}px`,
                    width: `${layout.labelWidth}px`,
                    height: `${layout.plotHeight}px`,
                    gridTemplateRows: `repeat(${visibleCategories.length}, minmax(0, 1fr))`,
                }}
            >
                {visibleCategories.map((category, rowIndex) => (
                    <div
                        key={`row-label-${category.categoryIndex}`}
                        className="flex items-center justify-end pr-3 text-right text-xs"
                        title={category.label}
                        style={{
                            zIndex:
                                hoveredCategory?.categoryIndex === category.categoryIndex ||
                                highlightedRowIndices.has(rowIndex)
                                    ? 1
                                    : 0,
                        }}
                    >
                        <span
                            className="max-w-full truncate rounded-sm px-1 py-0.5"
                            style={{
                                backgroundColor:
                                    hoveredCategory?.categoryIndex === category.categoryIndex
                                        ? "rgba(15, 23, 42, 0.78)"
                                        : highlightedRowIndices.has(rowIndex)
                                        ? "rgba(15, 23, 42, 0.78)"
                                        : "transparent",
                                color:
                                    hoveredCategory?.categoryIndex === category.categoryIndex ||
                                    highlightedRowIndices.has(rowIndex)
                                        ? "white"
                                        : undefined,
                                maxWidth:
                                    hoveredCategory?.categoryIndex === category.categoryIndex ? "none" : undefined,
                                overflow:
                                    hoveredCategory?.categoryIndex === category.categoryIndex ? "visible" : undefined,
                            }}
                        >
                            {category.label === "" ? "none" : category.label}
                        </span>
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
                {densityFields.map((field, columnIndex) => {
                    const color = getFieldColor(field.field);
                    return (
                        <div
                            key={`column-label-${field.field}`}
                            className="relative overflow-visible"
                            style={{
                                zIndex:
                                    hoveredField?.field === field.field || highlightedColumnIndices.has(columnIndex)
                                        ? 1
                                        : 0,
                            }}
                        >
                            <div
                                className="absolute bottom-2 left-2 whitespace-nowrap text-xs"
                                style={{
                                    color: `rgb(${color[0]}, ${color[1]}, ${color[2]})`,
                                    transform: "rotate(-35deg)",
                                    transformOrigin: "left bottom",
                                    backgroundColor:
                                        hoveredField?.field === field.field || highlightedColumnIndices.has(columnIndex)
                                            ? "rgba(15, 23, 42, 0.78)"
                                            : "transparent",
                                    borderRadius: "2px",
                                    padding:
                                        hoveredField?.field === field.field || highlightedColumnIndices.has(columnIndex)
                                            ? "2px 4px"
                                            : undefined,
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
                    className="absolute inset-0 pointer-events-none"
                >
                    <div
                        className="absolute inset-0"
                        style={{
                            border: `1px solid ${gridLineColor}`,
                        }}
                    />
                    {Array.from({ length: Math.max(densityFields.length - 1, 0) }, (_, index) => (
                        <div
                            key={`vline-${densityFields[index + 1]?.field ?? index}`}
                            className="absolute top-0 bottom-0"
                            style={{
                                left: `${(index + 1) * layout.cellWidth}px`,
                                width: "1px",
                                backgroundColor: gridLineColor,
                            }}
                        />
                    ))}
                    {Array.from({ length: Math.max(visibleCategories.length - 1, 0) }, (_, index) => (
                        <div
                            key={`hline-${visibleCategories[index + 1]?.categoryIndex ?? index}`}
                            className="absolute left-0 right-0"
                            style={{
                                top: `${(index + 1) * layout.cellHeight}px`,
                                height: "1px",
                                backgroundColor: gridLineColor,
                            }}
                        />
                    ))}
                </div>
                <div
                    className="absolute inset-0 grid"
                    style={{
                        gridTemplateColumns: `repeat(${densityFields.length}, minmax(0, 1fr))`,
                        gridTemplateRows: `repeat(${visibleCategories.length}, minmax(0, 1fr))`,
                    }}
                >
                    {visibleCategories.flatMap((category, rowIndex) =>
                        densityFields.map((field, columnIndex) => {
                            const highlighted = isCellHighlighted(rowIndex, columnIndex);
                            return (
                            <button
                                key={`grid-${category.categoryIndex}-${field.field}`}
                                type="button"
                                onMouseEnter={() => {
                                    setHoveredCell({ rowIndex, columnIndex });
                                }}
                                onMouseLeave={() => {
                                    setHoveredCell((current) =>
                                        current?.rowIndex === rowIndex && current?.columnIndex === columnIndex
                                            ? null
                                            : current,
                                    );
                                }}
                                onClick={(event) => {
                                    updateHighlightedCells(rowIndex, columnIndex, event.shiftKey);
                                }}
                                onKeyDown={(event) => {
                                    if (event.key === "Escape") {
                                        event.preventDefault();
                                        event.stopPropagation();
                                        clearHighlightedCells();
                                        return;
                                    }
                                    if (event.key !== "Enter" && event.key !== " ") return;
                                    event.preventDefault();
                                    updateHighlightedCells(rowIndex, columnIndex, event.shiftKey);
                                }}
                                aria-label={`${categoryColumn.name}: ${
                                    category.label === "" ? "none" : category.label
                                }, ${field.name === "" ? "none" : field.name}`}
                                aria-pressed={highlighted}
                                title={`${category.label === "" ? "none" : category.label} / ${field.name === "" ? "none" : field.name}`}
                                className="m-0 appearance-none border-0 p-0"
                                style={{
                                    backgroundColor:
                                        highlighted
                                            ? "rgba(255, 255, 255, 0.08)"
                                            : hoveredRowIndex === rowIndex && hoveredColumnIndex === columnIndex
                                            ? "rgba(255, 255, 255, 0.04)"
                                            : "transparent",
                                    boxShadow:
                                        highlighted
                                            ? "inset 0 0 0 1px rgba(255, 255, 255, 0.4)"
                                            : hoveredRowIndex === rowIndex && hoveredColumnIndex === columnIndex
                                            ? "inset 0 0 0 1px rgba(255, 255, 255, 0.22)"
                                            : "none",
                                    cursor: "pointer",
                                }}
                            />
                        );
                    }),
                    )}
                </div>
                {hoveredCategory && hoveredField ? (
                    <div
                        className="absolute left-3 top-3 pointer-events-none rounded-md px-3 py-2 text-xs"
                        style={{
                            backgroundColor: "rgba(15, 23, 42, 0.82)",
                            color: "white",
                            maxWidth: `min(${Math.max(layout.labelWidth - 24, 180)}px, calc(100% - 24px))`,
                            boxShadow: "0 4px 16px rgba(15, 23, 42, 0.18)",
                        }}
                    >
                        <div className="break-words">
                            <strong>{categoryColumn.name}:</strong>{" "}
                            {hoveredCategory.label === "" ? "none" : hoveredCategory.label}
                        </div>
                        <div className="break-words">
                            {hoveredField.name === "" ? "none" : hoveredField.name}
                        </div>
                    </div>
                ) : null}
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
