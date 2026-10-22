import type { DataColumn, DataType, LoadedDataColumn, NumberDataType } from "@/charts/charts";
import { canvasMeasureText, stripSharedSubgroupSuffix } from "@/charts/dotplot/xAxisLabels";
import { useChart, useDataStore } from "@/react/context";
import { useChartSize, useConfig, useFieldSpec, useFilteredIndices, useParamColumnsExperimental } from "@/react/hooks";
import { AxisLeft } from "@visx/axis";
import { scaleBand } from "@visx/scale";
import { observer } from "mobx-react-lite";
import { useEffect, useMemo, useState } from "react";
import DotPlotPanel, { type DotPlotHover, type DotPlotLabelHover } from "./DotPlotPanel";
import type { DotPlotReact } from "./DotPlotReactWrapper";
import { type DotPlotReactConfig, getColorRange } from "./dotPlotConfig";
import { DOT_PLOT_MAX_RADIUS, getDotPlotLayout, getVisibleRows, getXLabelLayout } from "./dotPlotLayout";
import { type LegendLayout, getLegendLayout } from "./dotPlotLegends";
import { MISSING_VALUE_COLOR, binnedColorScale, getDotPlotPalette, toCssColors } from "./dotPlotPalettes";
import { type DotPlotStats, queryDotPlotStats } from "./dotPlotStats";

type CategoryColumn = LoadedDataColumn<"text" | "text16">;
type NumberColumn = LoadedDataColumn<NumberDataType>;

function isCategoryColumn(column: DataColumn<DataType> | undefined): column is CategoryColumn {
    return column?.data !== undefined && (column.datatype === "text" || column.datatype === "text16");
}

function isNumberColumn(column: LoadedDataColumn<DataType>): column is NumberColumn {
    return column.datatype === "double" || column.datatype === "integer" || column.datatype === "int32";
}

/** Everything the stats depend on apart from the filtered rows. */
type StatsInputs = {
    rowColumn: CategoryColumn;
    fieldColumns: NumberColumn[];
    split?: { column: CategoryColumn; panels: number[] };
    threshold?: number;
};

type StatsState = {
    inputs: StatsInputs;
    stats: DotPlotStats | null;
    /** Colour range from the first stats for these inputs, kept while the filters change. */
    autoRange: [number, number] | null;
    error: string | null;
};

const NO_STATS = { stats: null, autoRange: null, error: null };

/**
 * Recomputes the stats in a worker whenever the inputs or filtered rows change. While a filter
 * change is being recomputed the previous stats are returned, but never stats for other inputs.
 */
function useDotPlotStats(
    inputs: StatsInputs | undefined,
    rows: Uint32Array,
): Pick<StatsState, "stats" | "autoRange" | "error"> {
    const [state, setState] = useState<StatsState | null>(null);
    useEffect(() => {
        if (!inputs) return;
        const { rowColumn, fieldColumns, split, threshold } = inputs;
        const controller = new AbortController();
        queryDotPlotStats(
            {
                rowCodes: rowColumn.data,
                rowCount: rowColumn.values.length,
                fields: fieldColumns.map((c) => c.data),
                rows,
                split: split && { codes: split.column.data, panels: split.panels },
                threshold,
            },
            controller.signal,
        )
            .then((stats) =>
                setState((prev) => ({
                    inputs,
                    stats,
                    autoRange: prev?.inputs === inputs && prev.autoRange ? prev.autoRange : stats.meanRange,
                    error: null,
                })),
            )
            .catch((error: unknown) => {
                if (controller.signal.aborted) return;
                const message = error instanceof Error ? error.message : String(error);
                setState({ inputs, stats: null, autoRange: null, error: message });
            });
        return () => controller.abort();
    }, [inputs, rows]);
    if (!state || state.inputs !== inputs) return NO_STATS;
    return state;
}

function Message({ children }: { children: string }) {
    return <div className="flex h-full items-center justify-center text-sm">{children}</div>;
}

const DotPlotComponent = observer(function DotPlotComponent() {
    const chart = useChart<DotPlotReactConfig, DotPlotReact>();
    const config = useConfig<DotPlotReactConfig>();
    const dataStore = useDataStore();
    const [width, height] = useChartSize();
    const columns = useParamColumnsExperimental();
    const rows = useFilteredIndices();
    const splitSpec = useFieldSpec(config.split_by ?? undefined);
    const [hovered, setHovered] = useState<DotPlotHover | null>(null);
    const [hoveredLabel, setHoveredLabel] = useState<DotPlotLabelHover | null>(null);

    // param is [rowColumn, ...fields]
    const { rowColumn, fieldColumns } = useMemo(() => {
        const [first, ...rest] = columns;
        return {
            rowColumn: isCategoryColumn(first) ? first : undefined,
            fieldColumns: rest.filter(isNumberColumn),
        };
    }, [columns]);

    // keyed by value: config.split_values is replaced, not mutated, but may be an equal copy
    const splitValuesKey = JSON.stringify(config.split_values);
    const splitColumn = config.split_by && isCategoryColumn(splitSpec) ? splitSpec : undefined;
    const panelCodes = useMemo(() => {
        if (!splitColumn) return null;
        const values: string[] = JSON.parse(splitValuesKey);
        return values.map((v) => splitColumn.values.indexOf(v)).filter((code) => code >= 0);
    }, [splitColumn, splitValuesKey]);

    const threshold = config.threshold_value;
    const inputs = useMemo<StatsInputs | undefined>(() => {
        if (!rowColumn) return undefined;
        if (config.split_by && (!splitColumn || !panelCodes?.length)) return undefined;
        const split = splitColumn && panelCodes ? { column: splitColumn, panels: panelCodes } : undefined;
        return { rowColumn, fieldColumns, split, threshold };
    }, [rowColumn, fieldColumns, config.split_by, splitColumn, panelCodes, threshold]);
    const { stats, autoRange, error } = useDotPlotStats(inputs, rows);

    useEffect(() => {
        chart.setAutoMeanRange(autoRange);
    }, [chart, autoRange]);

    const fieldLabels = useMemo(
        () =>
            stripSharedSubgroupSuffix(
                fieldColumns.map((c) => ({ field: c.field, name: c.name, inSubgroup: Boolean(c.subgroup) })),
            ),
        [fieldColumns],
    );

    const { palette, log, mode, min: manualMin, max: manualMax } = config.color_scale;
    const colorRange = getColorRange({ palette, log, mode, min: manualMin, max: manualMax }, autoRange);
    const rangeMin = colorRange?.[0];
    const rangeMax = colorRange?.[1];
    const { meanColor, legendColors } = useMemo(() => {
        if (rangeMin === undefined || rangeMax === undefined) {
            return { meanColor: () => MISSING_VALUE_COLOR, legendColors: null };
        }
        const colors = dataStore.getColumnColors(
            { datatype: "double" },
            {
                overideValues: {
                    min: rangeMin,
                    max: rangeMax,
                    colors: [...getDotPlotPalette(palette)],
                    colorLogScale: log,
                },
            },
        );
        return { meanColor: binnedColorScale(colors, rangeMin, rangeMax), legendColors: toCssColors(colors) };
    }, [dataStore, rangeMin, rangeMax, palette, log]);

    // Layout is worked out before the early returns so the legends can be updated from an effect.
    const visibleRows = stats && rowColumn ? getVisibleRows(stats.panels, rowColumn.values, config.y_axis_order) : [];
    const showColorLegend = config.color_legend.display;
    const showFractionLegend = config.fraction_legend.display;
    const panelTitles = splitColumn && panelCodes ? panelCodes.map((code) => splitColumn.values[code]) : null;
    const { x: xAxis, y: yAxis } = config.axis;
    const fracMax = stats?.fracMax;
    const sizes = {
        width,
        yAxisWidth: yAxis.size,
        titleHeight: panelTitles ? xAxis.tickfont + 12 : 0,
        panelGap: panelTitles ? 16 : 0,
        panelCount: stats?.panels.length ?? 0,
        fieldCount: fieldColumns.length,
        rowCount: visibleRows.length,
    };
    // Panel widths don't depend on the x labels, but the labels' height does depend on the widths.
    const across = getDotPlotLayout({ ...sizes, height, xAxisHeight: 0 });
    const plotLeft = across.panelLefts[0] ?? 0;
    const plotWidth = across.panelWidth;
    const { fit_labels, truncate_labels, thin_labels, rotate_labels, tickfont, size: xAxisSize } = xAxis;
    const xLabels = useMemo(
        () =>
            getXLabelLayout(
                { fit_labels, truncate_labels, thin_labels, rotate_labels, tickfont, size: xAxisSize },
                fieldLabels,
                { plotLeft, plotWidth, chartHeight: height },
                canvasMeasureText,
            ),
        [
            fit_labels,
            truncate_labels,
            thin_labels,
            rotate_labels,
            tickfont,
            xAxisSize,
            fieldLabels,
            plotLeft,
            plotWidth,
            height,
        ],
    );
    // The fraction legend's dots match the panels' dots, whose size depends on the space left
    // after the legend band: lay out for the largest dots, then again for the size that gives.
    const layoutBelow = (bandHeight: number) =>
        getDotPlotLayout({ ...sizes, height: height - bandHeight, xAxisHeight: xLabels.margin });
    const legendLayoutFor = (radius: number) =>
        getLegendLayout({
            width,
            height,
            colorRange:
                showColorLegend && rangeMin !== undefined && rangeMax !== undefined ? [rangeMin, rangeMax] : null,
            fracMax: showFractionLegend && fracMax !== undefined ? fracMax : null,
            maxRadius: radius,
        });
    const rawLegendLayout = legendLayoutFor(layoutBelow(legendLayoutFor(DOT_PLOT_MAX_RADIUS).height).maxRadius);
    const layout = layoutBelow(rawLegendLayout.height);
    const { maxRadius } = layout;

    // keyed by value, so the legends are only told about a real change
    const { height: bandHeight, colorWidth, fractionItems } = rawLegendLayout;
    const [colorLeft, colorTop] = rawLegendLayout.color ?? [];
    const [fractionLeft, fractionTop] = rawLegendLayout.fraction ?? [];
    const legendLayout = useMemo<LegendLayout>(
        () => ({
            height: bandHeight,
            colorWidth,
            fractionItems,
            color: colorLeft !== undefined && colorTop !== undefined ? [colorLeft, colorTop] : null,
            fraction: fractionLeft !== undefined && fractionTop !== undefined ? [fractionLeft, fractionTop] : null,
        }),
        [bandHeight, colorWidth, fractionItems, colorLeft, colorTop, fractionLeft, fractionTop],
    );
    useEffect(() => {
        chart.updateLegends(
            legendColors && rangeMin !== undefined && rangeMax !== undefined
                ? { colors: legendColors, range: [rangeMin, rangeMax], width: legendLayout.colorWidth }
                : null,
            fracMax !== undefined && maxRadius > 0 ? { fracMax, maxRadius, items: legendLayout.fractionItems } : null,
            legendLayout,
        );
    }, [chart, legendColors, rangeMin, rangeMax, fracMax, maxRadius, legendLayout]);

    if (error) return <Message>{`Couldn't compute the dot plot: ${error}`}</Message>;
    if (columns.length > 0 && !rowColumn) return <Message>Choose a category column for the y axis.</Message>;
    if (fieldColumns.length === 0 && rowColumn) return <Message>Choose one or more fields for the x axis.</Message>;
    if (config.split_by && splitColumn && !panelCodes?.length) {
        return <Message>Choose which values to split by.</Message>;
    }
    if (!stats || !rowColumn) return <Message>Loading…</Message>;
    if (visibleRows.length === 0) return <Message>No rows pass the current filters.</Message>;

    const yScale = scaleBand({ domain: visibleRows, range: [0, layout.plotHeight] });
    const hoveredCell = hovered && stats.panels[hovered.panel]?.[hovered.row]?.cells[hovered.field];

    return (
        <div className="relative h-full w-full">
            <svg
                width={width}
                height={height - bandHeight}
                className="absolute left-0"
                style={{ top: bandHeight }}
                role="img"
                aria-label={config.title}
            >
                <AxisLeft
                    left={layout.panelLefts[0]}
                    top={layout.top}
                    scale={yScale}
                    tickFormat={(row) => rowColumn.values[row]}
                    stroke="var(--text_color)"
                    tickStroke="var(--text_color)"
                    tickLabelProps={{ fill: "var(--text_color)", fontSize: yAxis.tickfont, textAnchor: "end" }}
                />
                {stats.panels.map((groups, panel) => (
                    <DotPlotPanel
                        // biome-ignore lint/suspicious/noArrayIndexKey: panels have a fixed order
                        key={panel}
                        panel={panel}
                        title={panelTitles?.[panel]}
                        groups={groups}
                        rows={visibleRows}
                        xLabels={xLabels}
                        layout={layout}
                        titleFontSize={tickfont + 2}
                        fracMax={stats.fracMax}
                        meanColor={meanColor}
                        hovered={hovered}
                        onHover={setHovered}
                        onLabelHover={setHoveredLabel}
                    />
                ))}
            </svg>
            {hovered && hoveredCell && (
                <div
                    className="pointer-events-none absolute rounded-md px-2 py-1 text-xs"
                    style={{
                        left: hovered.x + 12,
                        top: hovered.y + bandHeight + 12,
                        backgroundColor: "rgba(15, 23, 42, 0.85)",
                        color: "white",
                    }}
                >
                    <div>{fieldColumns[hovered.field]?.name}</div>
                    <div>
                        {rowColumn.name}: {rowColumn.values[hovered.row]}
                    </div>
                    {splitColumn && panelTitles && (
                        <div>
                            {splitColumn.name}: {panelTitles[hovered.panel]}
                        </div>
                    )}
                    <div>mean: {hoveredCell.mean.toFixed(2)}</div>
                    <div>expressing: {hoveredCell.frac.toFixed(0)}%</div>
                </div>
            )}
            {hoveredLabel && (
                // full name, as fitted labels may be shortened or left out
                <div
                    className="pointer-events-none absolute rounded-md px-2 py-1 text-xs"
                    style={{
                        left: hoveredLabel.x + 12,
                        top: hoveredLabel.y + bandHeight + 12,
                        backgroundColor: "rgba(15, 23, 42, 0.85)",
                        color: "white",
                    }}
                >
                    {fieldColumns[hoveredLabel.field]?.name}
                </div>
            )}
        </div>
    );
});

export default DotPlotComponent;
