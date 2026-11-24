import { action, makeObservable, observable } from "mobx";
import BaseChart, { type ColumnChangeEvent } from "@/charts/BaseChart";
import type { FieldName } from "@/charts/charts";
import getAxisGuiSpec from "@/charts/dialogs/utils/AxisSettingsGui";
import type DataStore from "@/datastore/DataStore";
import type Dimension from "@/datastore/Dimension";
import { type FieldSpec, flattenFields } from "@/lib/columnTypeHelpers";
import { g, isArray } from "@/lib/utils";
import { buildColorLegendSpec } from "@/react/legend/color_legend/buildColorLegendSpec";
import type { ColorLegendSpec } from "@/react/legend/color_legend/types";
import { BaseReactChart } from "../BaseReactChart";
import { buildChartSvg } from "../dotplot/dotPlotImage";
import MultiLineComponent from "./MultiLineComponent";
import {
    adaptMultiLineConfig,
    BAND_WIDTH_SLIDER,
    bandWidthFromFraction,
    getBandWidthFraction,
    getBrushFilterRange,
    getXDomain,
    MULTI_LINE_REACT_TYPE,
    type MultiLineConfig,
    type SavedMultiLineConfig,
    SCALE_TRIMS,
    type ScaleTrim,
} from "./multiLineConfig";

export class MultiLineReact extends BaseReactChart<MultiLineConfig> {
    /** Bumped when the data or colours of a column this chart draws are changed in place, so it redraws. */
    dataVersion = 0;
    /** Filters the other charts to the brushed range; this chart keeps drawing every row it would without it. */
    dim: Dimension;

    constructor(dataStore: DataStore, div: HTMLDivElement, originalConfig: SavedMultiLineConfig) {
        const config = adaptMultiLineConfig(originalConfig);
        const valueField = flattenFields(config.values)[0];
        if (!config.title && valueField) {
            config.title = dataStore.getColumnName(valueField) ?? valueField;
        }
        super(dataStore, div, config, MultiLineComponent);
        this.dim = this.dataStore.getDimension("range_dimension");
        makeObservable(this, {
            dataVersion: observable,
            bumpDataVersion: action,
            setBrush: action,
            setValueColumn: action,
            setCategories: action,
        });
    }

    /** Called by the brush; the component applies the filter, so a saved brush filters on load too. */
    setBrush(range: [number, number] | null) {
        this.config.brush = range;
    }

    /** Filters the other charts to the rows whose `field` is in the brushed range, or lifts the filter. */
    applyBrushFilter(field: FieldName | undefined, domain: [number, number] | undefined) {
        const { brush } = this.config;
        if (!brush || !field || !domain) {
            this.dim.removeFilter();
            this.resetButton.style.display = "none";
            return;
        }
        this.dim.filter("filterRange", [field], getBrushFilterRange(brush, domain));
        this.resetButton.style.display = "inline";
    }

    /** The chart menu's reset button. */
    removeFilter() {
        this.dim.removeFilter();
        this.setBrush(null);
    }

    onDataFiltered(dim?: unknown) {
        // "remove all filters" has already lifted this chart's filter
        if (dim === "all_removed" && this.config.brush) {
            this.setBrush(null);
            this.resetButton.style.display = "none";
        }
    }

    remove() {
        // also lifts the brush filter from the other charts
        this.dim.destroy();
        super.remove();
    }

    bumpDataVersion() {
        this.dataVersion++;
    }

    setValueColumn(field: FieldSpec) {
        const c = this.config;
        c.values = [field];
        // the brushed range was of the old column
        c.brush = null;
        // a saved label named the old column
        c.axis.x.label = null;
        // in data units: the old column's band width means nothing for the new one
        c.band_width = null;
    }

    setCategories(field: FieldName) {
        this.config.categories = field;
        this.setColorLegend();
    }

    /** For a link with `set_param`: `param` was `[value, categories]` in the legacy chart. */
    setParamAt(index: number, column: FieldName) {
        if (index === 0) this.setValueColumn(column);
        else if (index === 1) this.setCategories(column);
    }

    /** Every column the chart draws, as field names. */
    getFields(): FieldName[] {
        const { categories, values } = this.config;
        return [...(categories ? [categories] : []), ...flattenFields(values)];
    }

    // BaseChart only looks at `param`, which this chart leaves empty. Even when the edit refiltered,
    // the rows this chart draws may be unchanged, as they leave out its own brush
    onDataChanged(data: ColumnChangeEvent) {
        if (this.getFields().some((f) => data.columns.includes(f))) this.bumpDataVersion();
    }

    onColumnRemoved(column: FieldName) {
        if (!this.getFields().includes(column)) return false;
        this.remove();
        return true;
    }

    getColorLegendSpec(): ColorLegendSpec | null {
        const { categories } = this.config;
        if (!categories || !this.dataStore.columnIndex[categories]) return null;
        return buildColorLegendSpec(this.dataStore, categories);
    }

    /** For the chart menu's "create svg/png image". */
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

    /** The x domain of the value column for a trim (by default the current one), for the band width. */
    private getValueDomain(trim: ScaleTrim | null = this.config.scaletrim) {
        const field = flattenFields(this.config.values)[0];
        const column = field ? this.dataStore.columnIndex[field] : undefined;
        if (!column) return undefined;
        const quantiles = typeof column.quantiles === "object" ? column.quantiles : undefined;
        return getXDomain({ minMax: this.dataStore.getMinMaxForColumn(field), quantiles }, trim);
    }

    /** Column pickers bound to `categories` and `values`, in place of the ones `param` would give. */
    getParamSettings() {
        const c = this.config;
        return [
            g({
                type: "column",
                label: "Value (x axis)",
                current_value: c.values[0],
                columnType: "number",
                func: (v) => {
                    this.setValueColumn(isArray(v) ? v[0] : v);
                },
            }),
            g({
                type: "column",
                label: "Categories",
                current_value: c.categories ?? "",
                columnType: ["text", "text16"],
                func: (v) => {
                    const field = isArray(v) ? v[0] : v;
                    if (typeof field === "string") this.setCategories(field);
                },
            }),
        ];
    }

    getGeneralSettings() {
        const settings = super.getGeneralSettings();
        const params = settings.find((s) => s.type === "folder" && s.label === "Parameters");
        if (params?.type === "folder") params.current_value = this.getParamSettings();
        return settings;
    }

    /** The shared axis controls; their sizes apply when labels aren't fitted. */
    getAxisSettings() {
        const c = this.config;
        const folder = getAxisGuiSpec(c);
        folder.current_value.unshift(
            g({
                type: "check",
                label: "Fit margins to axis labels",
                current_value: c.fit_labels,
                func: (v) => {
                    c.fit_labels = v;
                },
            }),
        );
        return folder;
    }

    getSettings() {
        const c = this.config;
        const domain = this.getValueDomain();
        return [
            ...super.getSettings(),
            g({
                type: "check",
                label: "Stack Lines",
                current_value: c.stacked,
                func: (v) => {
                    c.stacked = v;
                },
            }),
            g({
                type: "check",
                label: "Fill Lines",
                current_value: c.fill,
                func: (v) => {
                    c.fill = v;
                },
            }),
            g({
                type: "radiobuttons",
                label: "Scale/Trim to Percentile",
                current_value: c.scaletrim ?? "none",
                choices: [["No Trim", "none"], ...SCALE_TRIMS.map((t): [string, string] => [t, t])],
                func: (v) => {
                    const trim = SCALE_TRIMS.find((t) => t === v) ?? null;
                    // band width is a fraction of the range when trimmed and in data units when not;
                    // a trim the column has no quantiles for draws untrimmed, so compare what's drawn
                    if (this.getValueDomain(trim)?.clamp !== this.getValueDomain()?.clamp) c.band_width = null;
                    c.scaletrim = trim;
                },
            }),
            ...(domain
                ? [
                      g({
                          // a fraction of the range either way, converted when set: the dialog keeps
                          // this control while the trim or value column changes under it
                          type: "slider",
                          label: "Band Width (fraction of range)",
                          ...BAND_WIDTH_SLIDER,
                          current_value: getBandWidthFraction(c, domain.domain, domain.clamp),
                          func: (v) => {
                              const now = this.getValueDomain();
                              if (now) c.band_width = bandWidthFromFraction(v, now.domain, now.clamp);
                          },
                      }),
                  ]
                : []),
            g({
                type: "slider",
                label: "Intervals",
                min: 10,
                max: 100,
                step: 1,
                current_value: c.intervals,
                func: (v) => {
                    c.intervals = v;
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
            this.getAxisSettings(),
        ];
    }
}

BaseChart.types[MULTI_LINE_REACT_TYPE] = {
    name: "Multi Line Chart",
    class: MultiLineReact,
    // only for the add-chart dialog's pickers: the chart moves them to `values` and `categories`
    params: [
        { type: "number", name: "Value (x axis)" },
        { type: ["text", "text16"], name: "Categories" },
    ],
};

export default 42;
