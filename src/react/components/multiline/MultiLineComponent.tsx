import { scaleLinear } from "d3-scale";
import { observer } from "mobx-react-lite";
import { useEffect, useMemo, useState } from "react";
import type { DataColumn, DataType, LoadedDataColumn, NumberDataType } from "@/charts/charts";
import { useChart, useDataStore } from "@/react/context";
import { useChartSize, useConfig, useFieldSpec, useFilteredIndices } from "@/react/hooks";
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

/** The x domain and clamp for a column and trim, from its metadata. */
function useXDomain(
    dataStore: ReturnType<typeof useDataStore>,
    column: NumberColumn | undefined,
    trim: ScaleTrim | null,
) {
    return useMemo(() => {
        if (!column) return undefined;
        const minMax = dataStore.getMinMaxForColumn(column.field);
        const quantiles = typeof column.quantiles === "object" ? column.quantiles : undefined;
        return getXDomain({ minMax, quantiles }, trim);
    }, [dataStore, column, trim]);
}

const MultiLineComponent = observer(function MultiLineComponent() {
    const chart = useChart<MultiLineConfig, MultiLineReact>();
    const config = useConfig<MultiLineConfig>();
    const dataStore = useDataStore();
    const [width, height] = useChartSize();
    const rows = useFilteredIndices();
    const categorySpec = useFieldSpec(config.categories ?? undefined);
    // one value for now; several, as a grid, come later
    const valueSpec = useFieldSpec(config.values[0]);
    const categoryColumn = isCategoryColumn(categorySpec) ? categorySpec : undefined;
    const valueColumn = isNumberColumn(valueSpec) ? valueSpec : undefined;

    const { scaletrim, band_width, intervals, stacked, fill, axis } = config;
    const xDomain = useXDomain(dataStore, valueColumn, scaletrim);
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
    }, [categoryColumn, valueColumn, xDomain, intervals, band_width]);
    const state = useDensities(inputs, rows, chart.dataVersion);

    const colors = useMemo(
        () => (categoryColumn ? dataStore.getColumnColors(categoryColumn.field).map(String) : []),
        [dataStore, categoryColumn],
    );

    // the legend shows the category column, which may only now have loaded or have changed
    // biome-ignore lint/correctness/useExhaustiveDependencies: runs when the category column changes
    useEffect(() => {
        chart.setColorLegend();
    }, [chart, categoryColumn]);

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
            axis={axis}
            xLabel={axis.x.label ?? valueColumn.name}
            yLabel={axis.y.label ?? "density"}
            onCurveHover={(event, category) => chart.showToolTip(event.nativeEvent, escapeHtml(labels[category] ?? ""))}
            onCurveLeave={() => chart.hideToolTip()}
        />
    );
});

export default MultiLineComponent;
