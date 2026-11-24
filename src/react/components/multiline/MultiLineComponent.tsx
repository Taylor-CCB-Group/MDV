import { scaleLinear } from "d3-scale";
import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useMemo, useState } from "react";
import type { DataColumn, DataType, LoadedDataColumn, NumberDataType } from "@/charts/charts";
import { useChart, useDataStore } from "@/react/context";
import { useChartSize, useConfig, useFieldSpec, useOwnedFilteredIndices } from "@/react/hooks";
import { escapeHtml } from "@/utilities/Utilities";
import MultiLinePanel from "./MultiLinePanel";
import type { MultiLineReact } from "./MultiLineReactWrapper";
import { getBandwidth, getXDomain, type MultiLineConfig, type ScaleTrim } from "./multiLineConfig";
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
 * The x domain and clamp for a column and trim, from its metadata. Kept by value, so it changes
 * only when the domain does: the densities are discarded rather than redrawn when it changes.
 */
function useXDomain(
    dataStore: ReturnType<typeof useDataStore>,
    column: NumberColumn | undefined,
    trim: ScaleTrim | null,
    dataVersion: number,
) {
    // biome-ignore lint/correctness/useExhaustiveDependencies: dataVersion changes when the column's data is replaced in place
    const xDomain = useMemo(() => {
        if (!column) return undefined;
        const minMax = dataStore.getMinMaxForColumn(column.field);
        const quantiles = typeof column.quantiles === "object" ? column.quantiles : undefined;
        return getXDomain({ minMax, quantiles }, trim);
    }, [dataStore, column, trim, dataVersion]);
    const lo = xDomain?.domain[0];
    const hi = xDomain?.domain[1];
    const clamp = xDomain?.clamp;
    return useMemo(
        () =>
            lo === undefined || hi === undefined || clamp === undefined
                ? undefined
                : { domain: [lo, hi] satisfies [number, number], clamp },
        [lo, hi, clamp],
    );
}

const MultiLineComponent = observer(function MultiLineComponent() {
    const chart = useChart<MultiLineConfig, MultiLineReact>();
    const config = useConfig<MultiLineConfig>();
    const dataStore = useDataStore();
    const [width, height] = useChartSize();
    // rows passing every filter but this chart's own brush, so its curves stay whole while brushing
    const { ownerVisibleRows: rows } = useOwnedFilteredIndices(chart.dim);
    const categorySpec = useFieldSpec(config.categories ?? undefined);
    // one value for now; several, as a grid, come later
    const valueSpec = useFieldSpec(config.values[0]);
    const categoryColumn = isCategoryColumn(categorySpec) ? categorySpec : undefined;
    const valueColumn = isNumberColumn(valueSpec) ? valueSpec : undefined;

    const { scaletrim, band_width, intervals, stacked, fill, axis } = config;
    // edits replace a column's data, and set its colours, on the same column object
    const { dataVersion } = chart;
    const xDomain = useXDomain(dataStore, valueColumn, scaletrim, dataVersion);
    // biome-ignore lint/correctness/useExhaustiveDependencies: dataVersion changes when a column's data is replaced in place
    const inputs = useMemo<DensityInputs | undefined>(() => {
        if (!categoryColumn || !valueColumn || !xDomain) return undefined;
        const { domain, clamp } = xDomain;
        return {
            categoryCodes: categoryColumn.data,
            categoryCount: categoryColumn.values.length,
            values: valueColumn.data,
            ticks: scaleLinear().domain(domain).ticks(intervals),
            bandwidth: getBandwidth({ band_width }, domain, clamp),
            clamp: clamp ? domain : undefined,
            domain,
        };
    }, [categoryColumn, valueColumn, xDomain, intervals, band_width, dataVersion]);
    const state = useDensities(inputs, rows, dataVersion);

    // the brush filters the other charts; applied from the config, so a saved brush filters on load
    const brushKey = JSON.stringify(config.brush);
    const valueField = valueColumn?.field;
    const domain = xDomain?.domain;
    // biome-ignore lint/correctness/useExhaustiveDependencies: config.brush is read by applyBrushFilter, keyed by value
    useEffect(() => {
        chart.applyBrushFilter(valueField, domain);
    }, [chart, brushKey, valueField, domain]);
    const setBrush = useCallback((range: [number, number] | null) => chart.setBrush(range), [chart]);

    // biome-ignore lint/correctness/useExhaustiveDependencies: dataVersion changes when the column's colours are set
    const colors = useMemo(
        () => (categoryColumn ? dataStore.getColumnColors(categoryColumn.field).map(String) : []),
        [dataStore, categoryColumn, dataVersion],
    );

    // the legend shows the category column, which may only now have loaded, changed or been recoloured
    // biome-ignore lint/correctness/useExhaustiveDependencies: runs when the category column or its colours change
    useEffect(() => {
        chart.setColorLegend();
    }, [chart, categoryColumn, dataVersion]);

    if (!config.categories || config.values.length === 0) return <Message>Choose categories and a value</Message>;
    if (!inputs) return null;
    if (state?.error) return <Message>{`Couldn't compute densities: ${state.error}`}</Message>;
    if (!state?.densities || !categoryColumn || !valueColumn) return null;
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
            xLabel={axis.x.label ?? valueColumn.name}
            yLabel={axis.y.label ?? "density"}
            onCurveHover={(event, category) => chart.showToolTip(event.nativeEvent, escapeHtml(labels[category] ?? ""))}
            onCurveLeave={() => chart.hideToolTip()}
            brush={{ value: config.brush, setValue: setBrush }}
        />
    );
});

export default MultiLineComponent;
