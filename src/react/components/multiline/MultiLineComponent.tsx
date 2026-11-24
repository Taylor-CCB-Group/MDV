import { scaleLinear } from "d3-scale";
import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useMemo, useState } from "react";
import type { DataColumn, DataType, LoadedDataColumn, NumberDataType } from "@/charts/charts";
import { useChart, useDataStore } from "@/react/context";
import {
    useChartID,
    useChartSize,
    useConfig,
    useFieldSpec,
    useFieldSpecs,
    useOwnedFilteredIndices,
} from "@/react/hooks";
import { escapeHtml } from "@/utilities/Utilities";
import { useChartArrayGrid } from "../../hooks/useChartArrayGrid";
import ChartArrayLayout from "../ChartArrayLayout";
import MultiLinePanel from "./MultiLinePanel";
import type { MultiLineReact } from "./MultiLineReactWrapper";
import {
    getBandwidth,
    getCellXDomains,
    type MultiLineConfig,
    type ScaleTrim,
    type ValueRange,
    type XDomain,
} from "./multiLineConfig";
import { type MultiLineDensities, type MultiLineDensityRequest, queryMultiLineDensities } from "./multiLineStats";

type CategoryColumn = LoadedDataColumn<"text" | "text16">;
type NumberColumn = LoadedDataColumn<NumberDataType>;

function isCategoryColumn(column: DataColumn<DataType> | undefined): column is CategoryColumn {
    return column?.data !== undefined && (column.datatype === "text" || column.datatype === "text16");
}

function isNumberColumn(column: DataColumn<DataType> | undefined): column is NumberColumn {
    return (
        column?.data !== undefined &&
        (column.datatype === "double" || column.datatype === "integer" || column.datatype === "int32")
    );
}

/** Everything the densities depend on apart from the filtered rows. */
type DensityInputs = Omit<MultiLineDensityRequest, "rows"> & { domain: [number, number] };

type DensityState = { inputs: DensityInputs; densities: MultiLineDensities | null; error: string | null };

function useDensities(inputs: DensityInputs | undefined, rows: Uint32Array, dataVersion: number) {
    const [state, setState] = useState<DensityState | null>(null);
    // biome-ignore lint/correctness/useExhaustiveDependencies: dataVersion changes when column data is edited in place
    useEffect(() => {
        if (!inputs) return;
        const controller = new AbortController();
        const { domain: _, ...request } = inputs;
        queryMultiLineDensities({ ...request, rows }, controller.signal)
            .then((densities) => setState({ inputs, densities, error: null }))
            .catch((error: unknown) => {
                if (controller.signal.aborted) return;
                const message = error instanceof Error ? error.message : String(error);
                setState({ inputs, densities: null, error: message });
            });
        return () => controller.abort();
    }, [inputs, rows, dataVersion]);
    // keep drawing the last result while the next one is computed, unless the inputs are different
    // in kind: another value column or trim (the old curves would be on the wrong axis) or another
    // category column (they would get the wrong labels and colours)
    if (!state || !inputs) return null;
    if (state.inputs.domain !== inputs.domain || state.inputs.categoryCodes !== inputs.categoryCodes) return null;
    return state;
}

function Message({ children }: { children: string }) {
    return <div className="flex h-full items-center justify-center text-sm">{children}</div>;
}

/**
 * The x domain of each value's cell, from column metadata (no data needs to be read). Kept by
 * value, so the domains change only when one does: a cell's densities are discarded rather than
 * redrawn when its domain changes.
 */
function useCellXDomains(
    dataStore: ReturnType<typeof useDataStore>,
    columns: NumberColumn[],
    trim: ScaleTrim | null,
    share: boolean,
    dataVersion: number,
) {
    // biome-ignore lint/correctness/useExhaustiveDependencies: dataVersion changes when a column's data is replaced in place
    const ranges = useMemo(
        () =>
            columns.map(
                (column): ValueRange => ({
                    minMax: dataStore.getMinMaxForColumn(column.field),
                    quantiles: typeof column.quantiles === "object" ? column.quantiles : undefined,
                }),
            ),
        [dataStore, columns, dataVersion],
    );
    const key = JSON.stringify([ranges, trim, share]);
    // biome-ignore lint/correctness/useExhaustiveDependencies: keyed by value
    return useMemo(() => getCellXDomains(ranges, trim, share), [key]);
}

type CellProps = {
    width: number;
    height: number;
    column: NumberColumn;
    categoryColumn: CategoryColumn;
    colors: readonly string[];
    xDomain: XDomain;
    /** Several cells: the axis label is the column's name, whatever the config's x label says. */
    inGrid: boolean;
};

/** One value's curves. Mounted only while its cell is in view, so only those cells compute densities. */
const MultiLineCell = observer(function MultiLineCell({
    width,
    height,
    column,
    categoryColumn,
    colors,
    xDomain,
    inGrid,
}: CellProps) {
    const chart = useChart<MultiLineConfig, MultiLineReact>();
    const config = useConfig<MultiLineConfig>();
    const { field } = column;
    // rows passing every filter but this value's own brush, so its curves stay whole while brushing;
    // a value never brushed has no dimension, and nothing of its own to leave out
    const { ownerVisibleRows: rows } = useOwnedFilteredIndices(chart.dims.get(field));
    const { band_width, intervals, stacked, fill, axis } = config;
    // edits replace a column's data on the same column object
    const { dataVersion } = chart;
    // biome-ignore lint/correctness/useExhaustiveDependencies: dataVersion changes when a column's data is replaced in place
    const inputs = useMemo<DensityInputs>(() => {
        const { domain, clamp } = xDomain;
        return {
            categoryCodes: categoryColumn.data,
            categoryCount: categoryColumn.values.length,
            values: column.data,
            ticks: scaleLinear().domain(domain).ticks(intervals),
            bandwidth: getBandwidth({ band_width }, domain, clamp),
            clamp: clamp ? domain : undefined,
            domain,
        };
    }, [categoryColumn, column, xDomain, intervals, band_width, dataVersion]);
    const state = useDensities(inputs, rows, dataVersion);
    const setBrush = useCallback((range: [number, number] | null) => chart.setBrush(field, range), [chart, field]);

    if (state?.error) return <Message>{`Couldn't compute densities: ${state.error}`}</Message>;
    if (!state?.densities) return null;
    const labels = categoryColumn.values;
    return (
        <MultiLinePanel
            width={width}
            height={height}
            domain={inputs.domain}
            ticks={state.inputs.ticks}
            curves={state.densities.curves}
            categoryLabels={labels}
            colors={colors}
            stacked={stacked}
            fill={fill}
            // copied, so this observer reads (and redraws on) the axis settings, which are edited in place
            axis={{ x: { ...axis.x }, y: { ...axis.y } }}
            fitLabels={config.fit_labels}
            xLabel={(inGrid ? null : axis.x.label) ?? column.name}
            yLabel={axis.y.label ?? "density"}
            onCurveHover={(event, category) => chart.showToolTip(event.nativeEvent, escapeHtml(labels[category] ?? ""))}
            onCurveLeave={() => chart.hideToolTip()}
            brush={{ value: config.brush[field] ?? null, setValue: setBrush }}
        />
    );
});

/** A cell's border, on each side (see chartArrayLayout.css). */
const CELL_BORDER = 1;

type GridProps = Omit<CellProps, "width" | "height" | "column" | "xDomain" | "inGrid"> & {
    columns: NumberColumn[];
    xDomains: XDomain[];
};

/** A cell per value, through `ChartArrayLayout`; cells out of view draw nothing. */
function MultiLineGrid({ columns, xDomains, ...cellProps }: GridProps) {
    const chartId = useChartID();
    const cells = useMemo(() => columns.map((c) => ({ key: c.field, label: c.name })), [columns]);
    const grid = useChartArrayGrid({ chartId, cells });
    const visible = useMemo(() => new Set(grid.visibleCellIndices), [grid.visibleCellIndices]);
    const renderCell = (index: number) => {
        const bounds = grid.getCellBounds(index);
        const column = columns[index];
        const xDomain = xDomains[index];
        if (!visible.has(index) || !bounds || !column || !xDomain) return null;
        return (
            <MultiLineCell
                {...cellProps}
                width={Math.max(0, Math.floor(bounds.width) - 2 * CELL_BORDER)}
                height={Math.max(0, Math.floor(bounds.height) - 2 * CELL_BORDER)}
                column={column}
                xDomain={xDomain}
                inGrid
            />
        );
    };
    return (
        <div ref={grid.scrollContainerRef} className="relative h-full w-full min-w-0 overflow-auto overflow-x-hidden">
            <ChartArrayLayout
                cellCount={grid.cellCount}
                cellKeys={grid.cellKeys}
                layoutRef={grid.setLayout}
                renderCell={renderCell}
            />
        </div>
    );
}

const MultiLineComponent = observer(function MultiLineComponent() {
    const chart = useChart<MultiLineConfig, MultiLineReact>();
    const config = useConfig<MultiLineConfig>();
    const dataStore = useDataStore();
    const [width, height] = useChartSize();
    const categorySpec = useFieldSpec(config.categories ?? undefined);
    const categoryColumn = isCategoryColumn(categorySpec) ? categorySpec : undefined;
    // every value's data, for now; loading only the cells in view is #693
    const valueSpecs = useFieldSpecs(config.values);
    const columns = useMemo(() => valueSpecs.filter(isNumberColumn), [valueSpecs]);

    // edits replace a column's data, and set its colours, on the same column object
    const { dataVersion } = chart;
    const xDomains = useCellXDomains(dataStore, columns, config.scaletrim, config.share_x, dataVersion);

    // the brushes filter the other charts; applied from the config, so saved brushes filter on load,
    // and for every value, not only those in view
    const brushKey = JSON.stringify(config.brush);
    const domainsKey = JSON.stringify(columns.map((c, i) => [c.field, xDomains[i]?.domain]));
    // biome-ignore lint/correctness/useExhaustiveDependencies: config.brush is read by applyBrushFilters; keyed by value
    useEffect(() => {
        chart.applyBrushFilters(new Map(columns.map((c, i) => [c.field, xDomains[i].domain])));
    }, [chart, brushKey, domainsKey]);

    // biome-ignore lint/correctness/useExhaustiveDependencies: dataVersion changes when the column's colours are set
    const colors = useMemo(
        () => (categoryColumn ? dataStore.getColumnColors(categoryColumn.field).map(String) : []),
        [dataStore, categoryColumn, dataVersion],
    );

    // the legend shows the category column, which may only now have loaded, changed or been
    // recoloured; one legend for the chart, as every cell has the same categories
    // biome-ignore lint/correctness/useExhaustiveDependencies: runs when the category column or its colours change
    useEffect(() => {
        chart.setColorLegend();
    }, [chart, categoryColumn, dataVersion]);

    if (!config.categories || config.values.length === 0) return <Message>Choose categories and values</Message>;
    if (!categoryColumn || columns.length === 0) return null;
    if (columns.length === 1) {
        return (
            <MultiLineCell
                width={width}
                height={height}
                column={columns[0]}
                categoryColumn={categoryColumn}
                colors={colors}
                xDomain={xDomains[0]}
                inGrid={false}
            />
        );
    }
    return <MultiLineGrid columns={columns} xDomains={xDomains} categoryColumn={categoryColumn} colors={colors} />;
});

export default MultiLineComponent;
