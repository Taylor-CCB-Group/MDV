import type { DataColumn, NumberDataType } from "@/charts/charts";
import type RangeDimension from "@/datastore/RangeDimension";
import * as d3 from "d3";
import { observer } from "mobx-react-lite";
import { useMemo, useEffect, useState, useCallback, useRef } from "react";
import { useDebounce } from "use-debounce";
import { useDimensionFilter, useConfig } from "../hooks";
import type { SelectionDialogConfig, RangeFilter } from "./SelectionDialogReact";

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

    const [histogram, setHistogram] = useState<number[]>([]);
    // this could be a more general utility function - expect to extract soon
    const queryHistogram = useCallback(async () => {
        const worker = new Worker(
            new URL("../../datastore/rawHistogramWorker.ts", import.meta.url),
        );
        worker.onmessage = (event) => {
            setHistogram(event.data);
            worker.terminate();
        };
        const isInt32 = column.datatype === "int32";
        const originalData = column.data as Float32Array | Int32Array;
        const data = new SharedArrayBuffer(originalData.length * 4);
        new Float32Array(data).set(originalData);
        worker.postMessage({
            data,
            min: column.minMax[0],
            max: column.minMax[1],
            bins: 100,
            isFloat: isInt32,
        });
    }, [column]);

    const [low, high] = value || column.minMax;
    const [min, max] = column.minMax;
    const lowFraction = (low - min) / (max - min);
    const highFraction = (high - min) / (max - min);

    return {
        value,
        step,
        histogram,
        lowFraction,
        highFraction,
        queryHistogram,
    };
}
// type set2d = ReturnType<typeof useState<[number, number]>>[1];
export type Range = [number, number];
export type set2d = (v: Range | null) => void; //nb, setting undefined can actually be problematic
export type FilterRangeType = ReturnType<typeof useRangeFilter>;
export type RangeProps = FilterRangeType & {
    setValue: set2d;
    minMax: Range;
    // probably want to review how these are specified / controlled
    histoWidth: number; //number of bins
    histoHeight: number; //height of the histogram
};

const useBrushX = (
    ref: React.RefObject<SVGSVGElement>,
    { value, setValue, minMax, histoWidth, histoHeight }: RangeProps, //consider different typing here
) => {
    const brushRef = useRef<ReturnType<typeof d3.brushX> | null>(null);
    // we need to be able to respond to changes in value - but without causing an infinite loop
    // or having the brush reset on every render
    const [initialValue] = useState(value);

    useEffect(() => {
        if (!ref.current) return;

        const svg = d3.select(ref.current);
        // Set up brush
        const brush = d3
            .brushX()
            .handleSize(1)
            .extent([
                [0, -2],
                [histoWidth, histoHeight + 2],
            ])
            .on("brush end", (event) => {
                if (event.selection) {
                    const [start, end] = event.selection.map((x: number) => {
                        if (!ref.current) {
                            console.error(
                                "No ref.current in brush event handler",
                            );
                            return 0;
                        }
                        const { width } = ref.current.getBoundingClientRect();
                        // Normalize x-coordinate to [minMax[0], minMax[1]]
                        const r = width / histoWidth;
                        const normalizedX = (r * x) / width;
                        return (
                            minMax[0] + normalizedX * (minMax[1] - minMax[0])
                        );
                    });
                    setValue([start, end]);
                } else {
                    // warning - the null value here does behave distinctly differently from undefined
                    // e.g. as of this writing, the reset button will be glitchy if we don't use null here
                    setValue(null); // null - reset to full range if brush is cleared
                }
            });

        brushRef.current = brush;

        // Apply the brush to the SVG
        const brushGroup = svg.append("g").attr("class", "brush").call(brush);
        // Initialize brush selection based on the initial value
        if (initialValue) {
            const [start, end] = initialValue.map(
                (v) => ((v - minMax[0]) / (minMax[1] - minMax[0])) * histoWidth,
            );
            brushGroup.call(brush.move, [start, end]); // Move the brush to the initial selection
        }

        // Apply `vectorEffect` directly to handles
        brushGroup
            .selectAll(".selection")
            .attr("vector-effect", "non-scaling-stroke");
        brushGroup
            .selectAll(".handle")
            .attr("vector-effect", "non-scaling-stroke");

        // Cleanup on unmount
        return () => {
            svg.select(".brush").remove();
        };
    }, [ref, setValue, minMax, histoWidth, histoHeight, initialValue]);

    const [debouncedValue] = useDebounce(value, 100, {
        equalityFn: (a, b) => {
            //although the type of input argument is [number, number] | null - they are undefined when component is unmounted
            //! which causes an exception here which breaks the whole chart
            //so rather than checking === null, we check for falsy values
            if (!a && !b) return true;
            if (!a || !b) return false;
            return a[0] === b[0] && a[1] === b[1];
        },
    });
    const setBrushValue = useCallback<set2d>(
        (v) => {
            if (!brushRef.current || !ref.current) return;
            const svg = d3.select(ref.current);

            if (!v) {
                // throw new Error("this is actually ok, but I want to test the error handling");
                //@ts-ignore life is too short
                svg.select(".brush").call(brushRef.current.move, null);
                return;
            }
            const [start, end] = v;
            const x0 =
                ((start - minMax[0]) / (minMax[1] - minMax[0])) * histoWidth;
            const x1 =
                ((end - minMax[0]) / (minMax[1] - minMax[0])) * histoWidth;
            //@ts-ignore life is too short
            svg.select(".brush").call(brushRef.current.move, [x0, x1]);
        },
        [minMax, histoWidth, ref],
    ); //why doesn't biome think we need brushRef?
    useEffect(() => {
        setBrushValue(debouncedValue);
    }, [debouncedValue, setBrushValue]);
};
export const Histogram = observer((props: RangeProps) => {
    const { histogram: data, queryHistogram, value } = props;
    const { histoWidth, histoHeight } = props;
    const ref = useRef<SVGSVGElement>(null);
    useBrushX(ref, props);
    const prefersDarkMode = window.mdv.chartManager.theme === "dark";
    const width = histoWidth;
    const height = histoHeight;
    const lineColor = prefersDarkMode ? "#fff" : "#000";
    // Find max value for vertical scaling
    const maxValue = Math.max(...data);

    // Define the padding and scaling factor
    const padding = 2;
    const xStep = data.length / (width + 1); // Space between points
    const yScale = (height - 2 * padding) / maxValue; // Scale based on max value

    const [hasQueried, setHasQueried] = useState(false);
    // if data changes, reset the hasQueried state
    useEffect(() => {
        data;
        setHasQueried(false);
    }, [data]);
    useEffect(() => {
        if (!ref.current) return;
        const observer = new IntersectionObserver(
            (entries) => {
                if (entries[0].isIntersecting && !hasQueried) {
                    setHasQueried(true);
                    queryHistogram();
                }
            },
            { rootMargin: "0px 0px 100px 0px" },
        );
        observer.observe(ref.current);
        // queryHistogram();
        return () => observer.disconnect();
    }, [queryHistogram, hasQueried]);

    // Generate the points for the polyline
    // ??? useMemo was wrong ????
    const points = useMemo(
        () =>
            data
                .map((value, index) => {
                    const x = index * xStep;
                    const y = height - padding - value * yScale;
                    return `${x},${y}`;
                })
                .join(" "),
        [data, xStep, yScale, height],
    );
    const v = value || props.minMax;
    return (
        <>
            <svg
                width={"100%"}
                height={height}
                viewBox={`0 0 ${width} ${height}`}
                preserveAspectRatio="none"
                ref={ref}
                cursor="move"
            >
                {/* Background polyline (the simple line connecting data points) */}
                <polyline
                    points={points}
                    fill="none"
                    stroke={lineColor}
                    strokeWidth="1.5"
                    // many thanks to ChatGPT for the following line (and the rest of the component
                    // but this would have been a real pain to figure out on my own)
                    vectorEffect="non-scaling-stroke" // Keeps the stroke width consistent
                />
                {/* d3.brushX will add more elements as a side-effect, handled in hook */}
            </svg>
            {/* <p className="flex justify-between"><em>{`${v[0].toFixed(2)}<`}</em> <em>{`<${v[1].toFixed(2)}`}</em></p> */}
        </>
    );
});
