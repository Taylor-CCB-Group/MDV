import { action, makeObservable, observable, runInAction } from "mobx";
import BaseChart, { type ColumnChangeEvent } from "@/charts/BaseChart";
import type { FieldName } from "@/charts/charts";
import getAxisGuiSpec from "@/charts/dialogs/utils/AxisSettingsGui";
import type DataStore from "@/datastore/DataStore";
import type Dimension from "@/datastore/Dimension";
import { type FieldSpecs, flattenFields } from "@/lib/columnTypeHelpers";
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
    getCellXDomains,
    keepBrushes,
    MULTI_LINE_REACT_TYPE,
    type MultiLineConfig,
    type SavedMultiLineConfig,
    SCALE_TRIMS,
    type ScaleTrim,
    withBrush,
} from "./multiLineConfig";

/**
 * An export draws the cells out of view this many at a time, as each computes its densities in a
 * worker of its own: a grid of hundreds of values would otherwise start hundreds at once.
 */
const CELL_EXPORT_BATCH = 16;
/** How long an export waits for a batch of cells; it moves on with what it has after. */
const CELL_EXPORT_TIMEOUT_MS = 20_000;

/** Resolves once the first `count` cells of `grid` have drawn something (curves, or a message), or on timeout. */
function waitForCells(grid: HTMLElement, count: number, timeoutMs: number) {
    const view = grid.ownerDocument.defaultView ?? window;
    const start = performance.now();
    return new Promise<void>((resolve) => {
        const check = () => {
            const cells = [...grid.querySelectorAll(".mdv-chart-array__cell")].slice(0, count);
            const drawn = cells.every((cell) => cell.childElementCount > 0);
            if (drawn || performance.now() - start > timeoutMs) resolve();
            else view.requestAnimationFrame(check);
        };
        // let the cells mount first
        view.requestAnimationFrame(check);
    });
}

export class MultiLineReact extends BaseReactChart<MultiLineConfig> {
    /** Bumped when the data or colours of a column this chart draws are changed in place, so it redraws. */
    dataVersion = 0;
    /**
     * A range dimension per brushed value, made on its first brush: each filters the other charts,
     * and the other cells, to the rows with its value in the brushed range. A cell draws every row
     * its own dimension excludes, so its curves stay whole while brushing.
     */
    readonly dims = observable.map<FieldName, Dimension>({}, { deep: false });
    /** The filter range each dimension has, as a string, so an unchanged brush isn't refiltered. */
    private appliedRanges = new Map<FieldName, string>();
    /**
     * While an image is made, a grid also draws its first this many cells, not only those in
     * view; 0 otherwise.
     */
    exportCellCount = 0;

    constructor(dataStore: DataStore, div: HTMLDivElement, originalConfig: SavedMultiLineConfig) {
        const config = adaptMultiLineConfig(originalConfig);
        const valueFields = flattenFields(config.values);
        if (!config.title && valueFields.length === 1) {
            config.title = dataStore.getColumnName(valueFields[0]) ?? valueFields[0];
        }
        super(dataStore, div, config, MultiLineComponent);
        makeObservable(this, {
            dataVersion: observable,
            bumpDataVersion: action,
            setBrush: action,
            removeFilter: action,
            exportCellCount: observable,
            setExportCellCount: action,
            setValues: action,
            setCategories: action,
        });
    }

    /** Called by a cell's brush; the component applies the filters, so saved brushes filter on load too. */
    setBrush(field: FieldName, range: [number, number] | null) {
        this.config.brush = withBrush(this.config.brush, field, range);
    }

    /**
     * Filters the other charts to each value's brushed range, or lifts the filter of a value with
     * no brush. `domains` are the x domains of the values whose data is loaded, by field.
     */
    applyBrushFilters(domains: ReadonlyMap<FieldName, [number, number]>) {
        const { brush } = this.config;
        const fields = flattenFields(this.config.values);
        for (const [field, dim] of [...this.dims]) {
            if (!fields.includes(field)) {
                // no longer drawn (a live query can drop it without `setValues`): a brush doesn't
                // outlive its column, or it would filter again if the column came back
                dim.destroy();
                runInAction(() => {
                    this.dims.delete(field);
                    this.config.brush = withBrush(this.config.brush, field, null);
                });
                this.appliedRanges.delete(field);
            } else if (!brush[field] && this.appliedRanges.has(field)) {
                dim.removeFilter();
                this.appliedRanges.delete(field);
            }
        }
        for (const field of fields) {
            const range = brush[field];
            const domain = domains.get(field);
            if (!range || !domain) continue;
            const filterRange = getBrushFilterRange(range, domain);
            const key = `${filterRange.min},${filterRange.max}`;
            if (this.appliedRanges.get(field) === key) continue;
            let dim = this.dims.get(field);
            if (!dim) {
                const created = this.dataStore.getDimension("range_dimension");
                runInAction(() => this.dims.set(field, created));
                dim = created;
            }
            dim.filter("filterRange", [field], filterRange);
            this.appliedRanges.set(field, key);
        }
        this.resetButton.style.display = this.appliedRanges.size ? "inline" : "none";
    }

    /** The chart menu's reset button: clears every brush. */
    removeFilter() {
        for (const dim of this.dims.values()) dim.removeFilter();
        this.appliedRanges.clear();
        this.config.brush = {};
        this.resetButton.style.display = "none";
    }

    onDataFiltered(dim?: unknown) {
        // "remove all filters" has already lifted every dimension's filter
        if (dim === "all_removed") {
            this.appliedRanges.clear();
            if (Object.keys(this.config.brush).length) runInAction(() => (this.config.brush = {}));
            this.resetButton.style.display = "none";
        }
    }

    remove() {
        // also lifts the brush filters from the other charts
        for (const dim of this.dims.values()) dim.destroy();
        super.remove();
    }

    bumpDataVersion() {
        this.dataVersion++;
    }

    setValues(values: FieldSpecs) {
        const c = this.config;
        const before = flattenFields(c.values);
        const after = flattenFields(values);
        c.values = values;
        // a brushed range only means something for its own column
        c.brush = keepBrushes(c.brush, after);
        if (before[0] !== after[0]) {
            // a saved label named the old column
            c.axis.x.label = null;
            // in data units: the old column's band width means nothing for the new one
            c.band_width = null;
        }
    }

    setCategories(field: FieldName) {
        this.config.categories = field;
        this.setColorLegend();
    }

    /** For a link with `set_param`: `param` was `[value, categories]` in the legacy chart. */
    setParamAt(index: number, column: FieldName) {
        if (index === 0) this.setValues([column]);
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
        void this.withEveryCellDrawn((container) => {
            const style = this.contentDiv.style;
            const textColor = style.getPropertyValue("--text_color");
            style.setProperty("--text_color", "black");
            const svg = buildChartSvg(container);
            if (textColor) style.setProperty("--text_color", textColor);
            else style.removeProperty("--text_color");
            if (type === "svg") {
                callback(new XMLSerializer().serializeToString(svg));
                return;
            }
            this.getImageFromSVG(svg, (canvas) => callback(canvas));
        });
    }

    setExportCellCount(count: number) {
        this.exportCellCount = count;
    }

    /**
     * Calls `draw` with the element to export. A grid only draws the cells in view, so for a grid
     * every cell is drawn first, a batch at a time (each computes its densities), and the grid's
     * whole layout is exported rather than the scrolled part of it. Once values load lazily (#693),
     * this also loads every value, which wants a limit or a warning for a large grid.
     */
    private async withEveryCellDrawn(draw: (container: HTMLElement) => void) {
        const grid = this.contentDiv.querySelector<HTMLElement>(".mdv-chart-array");
        if (!grid) {
            draw(this.contentDiv);
            return;
        }
        const cellCount = grid.querySelectorAll(".mdv-chart-array__cell").length;
        try {
            for (let count = CELL_EXPORT_BATCH; ; count += CELL_EXPORT_BATCH) {
                this.setExportCellCount(Math.min(count, cellCount));
                await waitForCells(grid, count, CELL_EXPORT_TIMEOUT_MS);
                if (count >= cellCount) break;
            }
            draw(grid);
        } finally {
            this.setExportCellCount(0);
        }
    }

    /**
     * The x domain of the first cell for a trim (by default the current one), for the band width.
     * Shared, it is every cell's.
     */
    private getValueDomain(trim: ScaleTrim | null = this.config.scaletrim) {
        const ranges = flattenFields(this.config.values).flatMap((field) => {
            const column = this.dataStore.columnIndex[field];
            if (!column) return [];
            const quantiles = typeof column.quantiles === "object" ? column.quantiles : undefined;
            return [{ minMax: this.dataStore.getMinMaxForColumn(field), quantiles }];
        });
        return getCellXDomains(ranges, trim, this.config.share_x)[0];
    }

    /** Column pickers bound to `categories` and `values`, in place of the ones `param` would give. */
    getParamSettings() {
        const c = this.config;
        return [
            g({
                type: "multicolumn",
                label: "Values (a cell each)",
                current_value: c.values,
                columnType: "number",
                func: (v) => {
                    this.setValues(isArray(v) ? [...v] : [v]);
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
        const bandWidthFraction = () => {
            const now = this.getValueDomain();
            return now ? getBandWidthFraction(c, now.domain, now.clamp) : BAND_WIDTH_SLIDER.min;
        };
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
            g({
                type: "check",
                label: "Share x axis between cells",
                current_value: c.share_x,
                func: (v) => {
                    c.share_x = v;
                },
            }),
            ...(domain
                ? [
                      g({
                          // the dialog keeps this control while the trim or value column changes
                          // under it, so it's a fraction of the range either way, converted when set,
                          // and it shows the config through a getter: the dialog's observable copy
                          // makes that a computed, so the slider follows a reset band width
                          type: "slider",
                          label: "Band Width (fraction of range)",
                          ...BAND_WIDTH_SLIDER,
                          get current_value() {
                              return bandWidthFraction();
                          },
                          set current_value(_) {
                              // `func` sets the config, which the getter reads
                          },
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
