import type { DataColumn, NumberDataType } from "@/charts/charts";
import type RangeDimension from "@/datastore/RangeDimension";
import { useMemo, useEffect, useState, useCallback, useRef, type SVGProps } from "react";
import { useDebounce } from "use-debounce";
import { useDimensionFilter, useConfig } from "../hooks";
import type { SelectionDialogConfig, RangeFilter } from "./SelectionDialogReact";
import { Brush } from "@visx/brush";
import { LinePath } from '@visx/shape';
import { ParentSize } from '@visx/responsive';
import { scaleLinear, scaleLog, scaleSymlog } from "@visx/scale";
import { observer } from "mobx-react-lite";
import type BaseBrush from "@visx/brush/lib/BaseBrush";
import { useHistogramQuery } from "../useHistogramQuery";
import { LinearProgress } from "@mui/material";


/**
 * This was exposed as a more general-purpose hook with useState,
 * moved to SelectionDialogComponent to handle state with mobx in the config.filters object.
 *
 * Now it's here, but it we still have some housekeeping to do.
 */
export function useRangeFilter(column: DataColumn<NumberDataType>) {
    const filter = useDimensionFilter(column) as RangeDimension;
    const filters = useConfig<SelectionDialogConfig>().filters;
    //nb - we may want to allow value to be null, rather than defaulting to minMax
    //relates to e.g. clearBrush function
    // const value = (filters[column.field] || column.minMax) as [number, number];
    const fVal = filters[column.field] as RangeFilter | null;
    // if (!isArray(fVal)) throw new Error("Expected range filter to be an array");
    const value = fVal;
    // const value = fVal;
    const isInteger = column.datatype.match(/int/);
    const { minMax } = column;
    const step = useMemo(() => {
        if (isInteger) return 1;
        // not sure this is totally correct - but there was a problem with very small ranges
        // this should be better...
        const small = Math.abs(minMax[1] - minMax[0]);
        return small < 0.001 ? small / 1000 : 0.001;
    }, [isInteger, minMax]);
    const [debouncedValue] = useDebounce(value, 10);

    // Effect to manage the filter state
    useEffect(() => {
        const value = debouncedValue;
        if (!value) {
            filter.removeFilter();
            return;
        }
        const [min, max] = value;
        filter.filter("filterRange", [column.field], { min, max }, true);
    }, [column, filter, debouncedValue]);

    const [low, high] = value || column.minMax;
    const [min, max] = column.minMax;
    const lowFraction = (low - min) / (max - min);
    const highFraction = (high - min) / (max - min);

    return {
        value,
        step,
        data: column.data || new Float32Array(),
        lowFraction,
        highFraction,
    };
}
// type set2d = ReturnType<typeof useState<[number, number]>>[1];
export type Range = [number, number];
export type set2d = (v: Range | null) => void; //nb, setting undefined can actually be problematic
export type FilterRangeType = ReturnType<typeof useRangeFilter>;
export type RangeProps = FilterRangeType & {
    setValue: set2d;
    domain: Range;
    // probably want to review how these are specified / controlled
    bins: number;
    histoHeight: number; //height of the histogram
    highlightValue?: number;
    xScaleType?: ScaleType;
    yScaleType?: ScaleType;
    name?: string;
    width?: number;
};
const ScaleTypes = {
    linear: scaleLinear,
    log: scaleLog,
    symlog: scaleSymlog,
} as const;
export type ScaleType = keyof typeof ScaleTypes;
const HistogramInner = observer(({
    xScaleType = "linear",
    yScaleType = "linear",
    width = 100,
    ...props
}: RangeProps) => {
    useEffect(() => {
        console.log("Histogram component mounted");
        return () => console.log("Histogram component unmounted");
    }, []);
    const { data: inputData, value } = props;
    const { bins, histoHeight, domain } = props;
    const ref = useRef<SVGSVGElement>(null);
    // nb, theme warrants use of observer
    const prefersDarkMode = window.mdv.chartManager.theme === "dark";
    const height = histoHeight;
    const lineColor = prefersDarkMode ? "#fff" : "#000";
    const selectedBrushStyle = useMemo<SVGProps<SVGRectElement>>(() => ({
        fill: prefersDarkMode ? "#fff" : "#000",
        stroke: prefersDarkMode ? "#fff" : "#000",
        strokeWidth: 1,
        fillOpacity: 0.1,
        vectorEffect: "non-scaling-stroke",
    }), [prefersDarkMode]);

    // todo lazy load
    // const [visible, setVisible] = useState(false);
    // useEffect(() => {
    //     if (!ref.current) return;
    //     const observer = new IntersectionObserver(
    //         (entries) => {
    //             if (entries[0].isIntersecting) {
    //                 setVisible(true);
    //             } else {
    //                 setVisible(false);
    //             }
    //         },
    //         { rootMargin: "0px 0px 100px 0px" },
    //     );
    //     observer.observe(ref.current);
    //     return () => observer.disconnect();
    // }, []);
    // Use TanStack Query for histogram data
    const {
        data = [],
        isLoading,
        isError,
        error
    } = useHistogramQuery(inputData, domain, bins, true, props.name);
    const maxValue = useMemo(() => Math.max(...data), [data]);
    // options for scales other than linear...
    // x log & y symlog seem to work for colour histograms...

    // todo use xScaleType for computing bins, rather than distorting the data
    // now they are linear 0 to data.length, mapping to [min, max]
    const histXScale = useMemo(() => scaleLinear({
        domain: [0, data.length - 1],
        range: domain,
    }), [domain, data.length]);
    const brushXScale = useMemo(() => ScaleTypes[xScaleType]({
        domain,
        range: [0, width],
    }), [domain, width, xScaleType]);
    const brushYScale = useMemo(() => ScaleTypes[yScaleType]({
        domain: [0, maxValue],
        range: [height, 0],
        nice: true,
    }), [height, maxValue, yScaleType]);
    const [initialValue] = useState(value);
    const initialBrushPosition = useMemo(() => {
        if (!initialValue) {
            return { start: { x: 0 }, end: { x: 0 } };
        }
        const start = { x: brushXScale(initialValue[0]) };
        const end = { x: brushXScale(initialValue[1]) };
        return { start, end };
    }, [brushXScale, initialValue]);
    // respond to reset button, which sets value to null.
    const brushRef = useRef<BaseBrush>(null);
    useEffect(() => {
        if (!brushRef.current) return;
        if (value === null) {
            brushRef.current.reset();
        }
    }, [value]);
    const highlightLine = useMemo(() => {
        if (props.highlightValue === undefined) return null;
        const x = brushXScale(props.highlightValue);
        const y = maxValue;
        return (
            <line
                x1={x}
                y1={0}
                x2={x}
                y2={y}
                stroke={lineColor}
                strokeWidth="1"
                opacity={0.5}
                strokeDasharray="5,5"
                vectorEffect="non-scaling-stroke"
            />
        );
    }, [brushXScale, maxValue, lineColor, props.highlightValue]);
    
    if (isLoading) return <div className="p-4"><LinearProgress /></div>;
    if (isError) return <div className="p-4">Error: {error?.message}</div>;

    return (
        <>
            <svg
                width={width}
                height={height}
                // with brushXScale.range = [0, width],
                // if I take viewBox out, it will be squished into left side,
                // but mouse events will be correct.
                // using ParentSize for responsiveness avoids this issue.
                // viewBox={`0 0 ${width} ${height}`}
                preserveAspectRatio="none"
                ref={ref}
                cursor="move"
            >
                {/* Background polyline (the simple line connecting data points) */}
                <LinePath
                    data={data}
                    x={(_, i) => brushXScale(histXScale(i))}
                    y={d => brushYScale(d)}
                    curve={undefined} // no curve, just a straight line
                    // points={points}
                    fill="none"
                    stroke={lineColor}
                    strokeWidth="1.5"
                    // many thanks to ChatGPT for the following line (and the rest of the component
                    // but this would have been a real pain to figure out on my own)
                    vectorEffect="non-scaling-stroke" // Keeps the stroke width consistent
                />
                <Brush
                    innerRef={brushRef}
                    xScale={brushXScale}
                    yScale={brushYScale}
                    width={width}
                    height={height}
                    brushDirection="horizontal"
                    initialBrushPosition={initialBrushPosition}
                    onChange={(v) => {
                        if (!v) return;
                        props.setValue([v.x0, v.x1]);
                    }}
                    useWindowMoveEvents //needs fixing wrt scale
                    selectedBoxStyle={selectedBrushStyle}
                />
                {highlightLine}
            </svg>
            {/* <p className="flex justify-between"><em>{`${v[0].toFixed(2)}<`}</em> <em>{`<${v[1].toFixed(2)}`}</em></p> */}
        </>
    );
});

export const Histogram = (props: RangeProps) => {
    const style = useMemo(() => ({
        width: "100%",
        height: props.histoHeight,
    }), [props.histoHeight]);
    return (
        <div style={style}>
            <ParentSize>
                {({ width, height }) => (
                    <HistogramInner
                        {...props}
                        // just need to sort out this business of what I mean by width
                        // bins={bins}
                        width={width}
                        histoHeight={height}
                    />
                )}
            </ParentSize>
        </div>
    )
}
