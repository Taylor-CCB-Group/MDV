import BaseChart from "@/charts/BaseChart";
import getAxisGuiSpec from "@/charts/dialogs/utils/AxisSettingsGui";
import type DataStore from "@/datastore/DataStore";
import type Dimension from "@/datastore/Dimension";
import { g } from "@/lib/utils";
import FractionLegend from "@/react/components/legend/FractionLegend";
import LegendWrapper from "@/react/components/legend/LegendWrapper";
import type { ColorLegendSpec } from "@/react/legend/color_legend/types";
import type { FractionLegendSpec } from "@/react/legend/fraction_legend/types";
import { makeObservable, observable, runInAction } from "mobx";
import { BaseReactChart } from "../BaseReactChart";
import DotPlotComponent from "./DotPlotComponent";
import type { DotPlotDot } from "./DotPlotPanel";
import {
    DOT_PLOT_REACT_OLD_TYPE,
    DOT_PLOT_REACT_TYPE,
    type DotPlotReactConfig,
    NO_COLUMN,
    adaptDotPlotConfig,
    getSplitFromDialog,
} from "./dotPlotConfig";
import { buildChartSvg } from "./dotPlotImage";
import { DOT_PLOT_LEGEND_MARGIN, type LegendLayout, buildFractionLegendSpec } from "./dotPlotLegends";

/** How many split values to show when a split column is first chosen: enough to compare two. */
const DEFAULT_SPLIT_VALUE_COUNT = 2;

const round2 = (x: number) => Math.round(x * 100) / 100;

/** `width` is the legend's content width, chosen to fit the chart (see `getLegendLayout`). */
export type DotPlotColorLegendInput = { colors: string[]; range: [number, number]; width: number };
/** `items` is the most dots to show, chosen to fit the chart (see `getLegendLayout`). */
export type DotPlotFractionLegendInput = { fracMax: number; maxRadius: number; items: number };

type FilterColumn = { field: string; data: ArrayLike<number> };
/** The columns a dot's cells are picked from: its rows, its field, and the split if there is one. */
export type DotPlotDotColumns = {
    rowColumn: FilterColumn;
    fieldColumn: FilterColumn;
    split?: { column: FilterColumn; code: number };
};

export class DotPlotReact extends BaseReactChart<DotPlotReactConfig> {
    /** The automatic colour range currently drawn, used to fill in the manual range. */
    autoMeanRange: [number, number] | null = null;
    fractionLegendWrapper = new LegendWrapper<FractionLegendSpec, DotPlotReactConfig>(this);
    private colorLegendInput: DotPlotColorLegendInput | null = null;
    private fractionLegendInput: DotPlotFractionLegendInput | null = null;
    /** Where the legends start, normally in the band the component keeps free above the panels. */
    private legendLayout: LegendLayout | null = null;
    /** Filters the other charts to a clicked dot's cells; this chart keeps showing all its dots. */
    dim: Dimension;
    /** The dot whose cells are filtered, if any. */
    selectedDot: DotPlotDot | null = null;

    constructor(dataStore: DataStore, div: HTMLDivElement, originalConfig: DotPlotReactConfig) {
        const config = adaptDotPlotConfig(originalConfig);
        const rowField = config.param[0];
        if (!config.title && typeof rowField === "string") {
            config.title = dataStore.getColumnName(rowField) ?? rowField;
        }
        super(dataStore, div, config, DotPlotComponent);
        this.dim = this.dataStore.getDimension("base_dimension");
        makeObservable(this, { selectedDot: observable.ref });
    }

    /**
     * Filters to a dot's cells, as the legacy dot plot does (`filterCatCol`): cells in its row
     * category (and split value) whose value for its field is above the threshold. Clicking the
     * filtered dot again clears the filter.
     */
    filterToDot(dot: DotPlotDot, { rowColumn, fieldColumn, split }: DotPlotDotColumns) {
        const current = this.selectedDot;
        if (current?.panel === dot.panel && current.row === dot.row && current.field === dot.field) {
            this.removeFilter();
            this.resetButton.style.display = "none";
            return;
        }
        const threshold = this.config.threshold_value ?? 0;
        const rows = rowColumn.data;
        const values = fieldColumn.data;
        const splitCodes = split?.column.data;
        const splitCode = split?.code;
        // false for NaN values, as in the stats
        const predicate = (i: number) =>
            rows[i] === dot.row && (splitCodes === undefined || splitCodes[i] === splitCode) && values[i] > threshold;
        const columns = [rowColumn.field, fieldColumn.field, ...(split ? [split.column.field] : [])];
        this.dim.filter("filterPredicate", columns, { predicate });
        runInAction(() => {
            this.selectedDot = dot;
        });
        this.resetButton.style.display = "inline";
    }

    /** Clears the dot filter, e.g. when the rows, fields or split change and the dot means something else. */
    clearDotFilter() {
        if (!this.selectedDot) return;
        this.removeFilter();
        this.resetButton.style.display = "none";
    }

    removeFilter() {
        this.dim.removeFilter();
        runInAction(() => {
            this.selectedDot = null;
        });
    }

    onDataFiltered(dim?: unknown) {
        // "remove all filters" has already cleared the dimension's filter
        if (dim === "all_removed") {
            runInAction(() => {
                this.selectedDot = null;
            });
            this.resetButton.style.display = "none";
        }
    }

    setAutoMeanRange(range: [number, number] | null) {
        this.autoMeanRange = range;
    }

    /**
     * Called by the component with what the legends should show and where they start. The legends
     * live outside the component's React tree (in `LegendWrapper` portals), so they're only
     * re-rendered on change. A null input keeps the current legend, e.g. while stats recompute.
     */
    updateLegends(
        color: DotPlotColorLegendInput | null,
        fraction: DotPlotFractionLegendInput | null,
        legendLayout: LegendLayout,
    ) {
        this.legendLayout = legendLayout;
        const prevColor = this.colorLegendInput;
        if (
            color &&
            (color.colors !== prevColor?.colors ||
                color.range[0] !== prevColor.range[0] ||
                color.range[1] !== prevColor.range[1] ||
                color.width !== prevColor.width)
        ) {
            this.colorLegendInput = color;
            this.setColorLegend();
        }
        const prevFraction = this.fractionLegendInput;
        if (
            fraction &&
            (fraction.fracMax !== prevFraction?.fracMax ||
                fraction.maxRadius !== prevFraction.maxRadius ||
                fraction.items !== prevFraction.items)
        ) {
            this.fractionLegendInput = fraction;
            this.showFractionLegend();
        }
    }

    /**
     * Removes both legends, e.g. while the split has no values to show; unlike a null input to
     * `updateLegends`, which keeps them while stats recompute. Their positions are kept.
     */
    clearLegends() {
        const legend = this.legend;
        if (legend) {
            runInAction(() => {
                this.config.color_legend.pos = [legend.offsetLeft, legend.offsetTop];
            });
            this.colorLegendWrapper.unmount();
            this.legend = undefined;
        }
        this.colorLegendInput = null;
        if (this.fractionLegendInput) {
            this.fractionLegendInput = null;
            runInAction(() => this.showFractionLegend());
        }
    }

    /** Brushing would filter cells by one field's values, which isn't what mean expression shows. */
    isColorLegendFilterable() {
        return false;
    }

    isColorLegendResizable() {
        return true;
    }

    getColorLegendSpec(): ColorLegendSpec | null {
        const input = this.colorLegendInput;
        if (!input) return null;
        const { colors, range, width } = input;
        return { kind: "continuous", label: "Mean expression", column: "", colors, range, width };
    }

    setColorLegend() {
        const cl = this.config.color_legend;
        if (cl.display) {
            // BaseChart calls this while the chart is still being built, before there's anything
            // to show; it would then save a default position over the panels
            if (!this.colorLegendInput) return;
            if (!cl.pos && !this.legend && this.legendLayout?.color) {
                cl.pos = [...this.legendLayout.color];
            }
        }
        super.setColorLegend();
    }

    /** The legacy dot plot's fraction legend, laid out in a row to fit the band above the panels. */
    showFractionLegend() {
        const c = this.config.fraction_legend;
        const current = this.fractionLegendWrapper.getWrapperElement();
        if (current) {
            c.position = [current.offsetLeft, current.offsetTop];
        }
        const input = this.fractionLegendInput;
        if (!c.display || !input) {
            this.fractionLegendWrapper.unmount();
            return;
        }
        const position = c.position ?? this.legendLayout?.fraction ?? [DOT_PLOT_LEGEND_MARGIN, DOT_PLOT_LEGEND_MARGIN];
        c.position = [...position];
        const spec = buildFractionLegendSpec(input.fracMax, input.maxRadius, input.items);
        this.fractionLegendWrapper.render(spec, { left: position[0], top: position[1] }, FractionLegend, {
            dragHandle: ".legend-drag-handle",
            resizable: true,
        });
    }

    /** For the chart menu's "create svg/png image": the plot with its legends, black text. */
    getImage(callback: (image: string | HTMLCanvasElement) => void, type: "png" | "svg") {
        const style = this.contentDiv.style;
        const textColor = style.getPropertyValue("--text_color");
        style.setProperty("--text_color", "black");
        const svg = buildChartSvg(this.contentDiv);
        if (textColor) style.setProperty("--text_color", textColor);
        else style.removeProperty("--text_color");

        if (type === "svg") {
            callback(new XMLSerializer().serializeToString(svg));
            return;
        }
        this.getImageFromSVG(svg, (canvas) => callback(canvas));
    }

    getConfig() {
        const el = this.fractionLegendWrapper.getWrapperElement();
        if (el) {
            this.config.fraction_legend.position = [el.offsetLeft, el.offsetTop];
        }
        return super.getConfig();
    }

    remove() {
        this.fractionLegendWrapper.unmount();
        // also lifts the dot filter from the other charts
        this.dim.destroy();
        super.remove();
    }

    getSplitSettings() {
        const c = this.config;
        const columns = this.dataStore.getColumnList(["text", "text16"], true);
        return [
            g({
                type: "dropdown",
                label: "Split by",
                current_value: c.split_by ?? "__none__",
                values: [columns, "name", "field"],
                func: (field) => {
                    const values: string[] = this.dataStore.columnIndex[field]?.values ?? [];
                    c.split_by = field === "__none__" ? null : field;
                    c.split_values = c.split_by ? values.slice(0, DEFAULT_SPLIT_VALUE_COUNT) : [];
                },
            }),
            g({
                type: "category_selection",
                label: "Split values",
                current_value: c.split_values,
                sourceColumn: () => c.split_by ?? undefined,
                getCurrentValue: () => c.split_values,
                func: (values) => {
                    c.split_values = values;
                },
            }),
        ];
    }

    getScaleSettings() {
        const cs = this.config.color_scale;
        const auto = this.autoMeanRange;
        return g({
            type: "folder",
            label: "Scale",
            current_value: [
                g({
                    type: "radiobuttons",
                    label: "Colour range",
                    current_value: cs.mode,
                    choices: [
                        ["Auto", "auto"],
                        ["Manual", "manual"],
                    ],
                    func: (v) => {
                        if (v !== "auto" && v !== "manual") return;
                        cs.mode = v;
                        // start the manual range from what's drawn now
                        if (v === "manual" && auto) {
                            cs.min ??= round2(auto[0]);
                            cs.max ??= round2(auto[1]);
                        }
                    },
                }),
                g({
                    type: "spinner",
                    label: "Manual min",
                    current_value: cs.min ?? (auto ? round2(auto[0]) : 0),
                    min: -1e9,
                    max: 1e9,
                    step: 0.1,
                    func: (v) => {
                        cs.min = v;
                    },
                }),
                g({
                    type: "spinner",
                    label: "Manual max",
                    current_value: cs.max ?? (auto ? round2(auto[1]) : 1),
                    min: -1e9,
                    max: 1e9,
                    step: 0.1,
                    func: (v) => {
                        cs.max = v;
                    },
                }),
                g({
                    type: "radiobuttons",
                    label: "Color Palette",
                    current_value: cs.palette,
                    choices: [
                        ["Reds", "reds"],
                        ["Blue-White-Red", "blue_white_red"],
                    ],
                    func: (v) => {
                        if (v === "reds" || v === "blue_white_red") cs.palette = v;
                    },
                }),
                g({
                    type: "check",
                    label: "Log Color Scale",
                    current_value: cs.log,
                    func: (v) => {
                        cs.log = v;
                    },
                }),
            ],
        });
    }

    getSettings() {
        const c = this.config;
        const settings = super.getSettings();
        // split by belongs with the other columns the chart is built from
        const params = settings.find((s) => s.type === "folder" && s.label === "Parameters");
        if (params?.type === "folder") {
            params.current_value.push(...this.getSplitSettings());
        }
        return [
            ...settings,
            g({
                type: "radiobuttons",
                label: "Y-axis order",
                current_value: c.y_axis_order,
                choices: [
                    ["Data order", "data"],
                    ["Alphabetical", "alphabetical"],
                ],
                func: (v) => {
                    if (v === "data" || v === "alphabetical") c.y_axis_order = v;
                },
            }),
            g({
                type: "check",
                label: "Show Color Legend",
                current_value: c.color_legend.display,
                func: (v) => {
                    c.color_legend.display = v;
                    this.setColorLegend();
                },
            }),
            g({
                type: "check",
                label: "Show Fraction Legend",
                current_value: c.fraction_legend.display,
                func: (v) => {
                    c.fraction_legend.display = v;
                    this.showFractionLegend();
                },
            }),
            g({
                type: "check",
                label: "Keep space above the plot for legends",
                current_value: c.reserve_legend_space,
                func: (v) => {
                    c.reserve_legend_space = v;
                },
            }),
            this.getScaleSettings(),
            this.getAxisSettings(),
        ];
    }

    /** The shared axis controls, after the legacy dot plot's x label options. */
    getAxisSettings() {
        const x = this.config.axis.x;
        const folder = getAxisGuiSpec(this.config, undefined, { rotateXLabels: true });
        folder.current_value.unshift(
            g({
                type: "check",
                label: "Auto-fit X axis labels",
                current_value: x.fit_labels,
                func: (v) => {
                    x.fit_labels = v;
                },
            }),
            g({
                type: "check",
                label: "Shorten long X axis labels",
                current_value: x.truncate_labels,
                func: (v) => {
                    x.truncate_labels = v;
                },
            }),
            g({
                type: "check",
                label: "Hide overlapping X axis labels",
                current_value: x.thin_labels,
                func: (v) => {
                    x.thin_labels = v;
                },
            }),
        );
        return folder;
    }
}

BaseChart.types[DOT_PLOT_REACT_TYPE] = {
    name: "Dot Plot (Compound)",
    class: DotPlotReact,
    params: [
        {
            type: ["text", "text16"],
            name: "Categories on y-axis",
        },
        {
            type: "_multi_column:number",
            name: "Fields on x axis",
        },
    ],
    extra_controls: (dataStore: DataStore) => [
        {
            type: "dropdown",
            name: "split_by",
            label: "Split by",
            values: [
                { name: "None", value: NO_COLUMN },
                ...dataStore
                    .getColumnList(["text", "text16"])
                    .map((c) => ({ name: c.name ?? "", value: c.field ?? "" })),
            ],
            defaultVal: NO_COLUMN,
        },
        {
            // starts empty: the values to compare are chosen, not guessed
            type: "category_selection",
            name: "split_values",
            label: "Split values",
            columnControl: "split_by",
            defaultVal: [],
        },
    ],
    // the dialog has already copied both controls onto the config; tidy them up
    init: (config: DotPlotReactConfig, dataStore: DataStore, extra: Record<string, unknown>) => {
        const field = typeof extra.split_by === "string" ? extra.split_by : undefined;
        const columnValues = field ? dataStore.columnIndex[field]?.values : undefined;
        Object.assign(config, getSplitFromDialog(extra.split_by, extra.split_values, columnValues));
    },
};

// Same chart under its earlier name, so views saved with it still open. Hidden from Add Chart,
// and named differently because the add-chart dialog looks chart types up by name.
BaseChart.types[DOT_PLOT_REACT_OLD_TYPE] = {
    ...BaseChart.types[DOT_PLOT_REACT_TYPE],
    name: "Dot Plot (Compound)",
    allow_user_add: false,
};

export default 42;
