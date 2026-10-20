import BaseChart from "@/charts/BaseChart";
import getAxisGuiSpec from "@/charts/dialogs/utils/AxisSettingsGui";
import type DataStore from "@/datastore/DataStore";
import { g } from "@/lib/utils";
import FractionLegend from "@/react/components/legend/FractionLegend";
import LegendWrapper from "@/react/components/legend/LegendWrapper";
import type { ColorLegendSpec } from "@/react/legend/color_legend/types";
import type { FractionLegendSpec } from "@/react/legend/fraction_legend/types";
import { BaseReactChart } from "../BaseReactChart";
import DotPlotComponent from "./DotPlotComponent";
import {
    DOT_PLOT_REACT_OLD_TYPE,
    DOT_PLOT_REACT_TYPE,
    type DotPlotReactConfig,
    NO_COLUMN,
    adaptDotPlotConfig,
    getSplitFromDialog,
} from "./dotPlotConfig";
import { DOT_PLOT_LEGEND_MARGIN, type LegendBand, buildFractionLegendSpec } from "./dotPlotLegends";

/** How many split values to show when a split column is first chosen: enough to compare two. */
const DEFAULT_SPLIT_VALUE_COUNT = 2;

const round2 = (x: number) => Math.round(x * 100) / 100;

export type DotPlotColorLegendInput = { colors: string[]; range: [number, number] };
export type DotPlotFractionLegendInput = { fracMax: number; maxRadius: number };

export class DotPlotReact extends BaseReactChart<DotPlotReactConfig> {
    /** The automatic colour range currently drawn, used to fill in the manual range. */
    autoMeanRange: [number, number] | null = null;
    fractionLegendWrapper = new LegendWrapper<FractionLegendSpec, DotPlotReactConfig>(this);
    private colorLegendInput: DotPlotColorLegendInput | null = null;
    private fractionLegendInput: DotPlotFractionLegendInput | null = null;
    /** Where the legends start, in the band the component keeps free above the panels. */
    private legendBand: LegendBand | null = null;

    constructor(dataStore: DataStore, div: HTMLDivElement, originalConfig: DotPlotReactConfig) {
        const config = adaptDotPlotConfig(originalConfig);
        const rowField = config.param[0];
        if (!config.title && typeof rowField === "string") {
            config.title = dataStore.getColumnName(rowField) ?? rowField;
        }
        super(dataStore, div, config, DotPlotComponent);
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
        band: LegendBand,
    ) {
        this.legendBand = band;
        const prevColor = this.colorLegendInput;
        if (
            color &&
            (color.colors !== prevColor?.colors ||
                color.range[0] !== prevColor.range[0] ||
                color.range[1] !== prevColor.range[1])
        ) {
            this.colorLegendInput = color;
            this.setColorLegend();
        }
        const prevFraction = this.fractionLegendInput;
        if (fraction && (fraction.fracMax !== prevFraction?.fracMax || fraction.maxRadius !== prevFraction.maxRadius)) {
            this.fractionLegendInput = fraction;
            this.showFractionLegend();
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
        return { kind: "continuous", label: "Mean expression", column: "", colors: input.colors, range: input.range };
    }

    setColorLegend() {
        const cl = this.config.color_legend;
        if (cl.display) {
            // BaseChart calls this while the chart is still being built, before there's anything
            // to show; it would then save a default position over the panels
            if (!this.colorLegendInput) return;
            if (!cl.pos && !this.legend && this.legendBand?.color) {
                cl.pos = [...this.legendBand.color];
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
        const position = c.position ?? this.legendBand?.fraction ?? [DOT_PLOT_LEGEND_MARGIN, DOT_PLOT_LEGEND_MARGIN];
        c.position = [...position];
        const spec = buildFractionLegendSpec(input.fracMax, input.maxRadius);
        this.fractionLegendWrapper.render(spec, { left: position[0], top: position[1] }, FractionLegend, {
            dragHandle: ".legend-drag-handle",
            resizable: true,
        });
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
        if (params && Array.isArray(params.current_value)) {
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
            this.getScaleSettings(),
            getAxisGuiSpec(c, undefined, { rotateXLabels: true }),
        ];
    }
}

BaseChart.types[DOT_PLOT_REACT_TYPE] = {
    name: "Compound Dot Plot",
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
            sourceControl: "split_by",
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
    name: "Compound Dot Plot (old type name)",
    allow_user_add: false,
};

export default 42;
