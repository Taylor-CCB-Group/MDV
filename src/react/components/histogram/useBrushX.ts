import * as d3 from "d3";
import { useCallback, useEffect, useRef, type RefObject } from "react";

export type Range = [number, number];
export type BrushXScaleType = "linear" | "log";

export type BrushXConfig = {
    value: Range | null;
    setValue: (value: Range | null) => void;
    minMax: Range;
};

export type BrushXLayout = {
    x?: number;
    y?: number;
    width: number;
    height: number;
    handleSize?: number;
};

export type BrushXStyle = {
    selectionFill?: string;
    selectionStroke?: string;
    handleFill?: string;
};

export type BrushXOptions = {
    layout: BrushXLayout;
    style?: BrushXStyle;
    updateOn?: "brush" | "end";
};

export const createBrushScale = (
    scaleType: BrushXScaleType,
    domain: Range,
    range: Range,
) =>
    scaleType === "log"
        ? d3.scaleSymlog().domain(domain).range(range)
        : d3.scaleLinear().domain(domain).range(range);

export function useBrushX(
    ref: RefObject<SVGSVGElement | null>,
    brushConfig: BrushXConfig | undefined,
    xScaleType: BrushXScaleType,
    options: BrushXOptions,
) {
    const brushRef = useRef<ReturnType<typeof d3.brushX> | null>(null);
    const isUserBrushingRef = useRef(false);
    const minMax = brushConfig?.minMax;
    const setValue = brushConfig?.setValue;
    const {
        x = 0,
        y = -2,
        width,
        height,
        handleSize = 1.25,
    } = options.layout;
    const {
        selectionFill = "rgba(15, 23, 42, 0.08)",
        selectionStroke = "rgba(15, 23, 42, 0.7)",
        handleFill = "rgba(15, 23, 42, 0.88)",
    } = options.style ?? {};
    const updateOn = options.updateOn ?? "brush";

    useEffect(() => {
        if (!ref.current || !minMax || !setValue) return;

        const svg = d3.select(ref.current);
        const xScale = createBrushScale(xScaleType, minMax, [x, x + width]);
        const brush = d3.brushX()
            .handleSize(handleSize)
            .extent([[x, y], [x + width, y + height]])
            .on("brush end", (event) => {
                if (!event.sourceEvent) return;
                if (event.type === "brush") {
                    isUserBrushingRef.current = true;
                } else if (event.type === "end") {
                    requestAnimationFrame(() => {
                        isUserBrushingRef.current = false;
                    });
                }
                if (updateOn === "end" && event.type !== "end") return;
                if (event.selection) {
                    const [start, end] = event.selection.map((position: number) =>
                        xScale.invert(position),
                    );
                    setValue([start, end]);
                } else {
                    setValue(null);
                }
            });

        brushRef.current = brush;
        const brushGroup = svg.append("g").attr("class", "brush").call(brush);
        brushGroup
            .selectAll(".selection")
            .attr("fill", selectionFill)
            .attr("stroke", selectionStroke)
            .attr("vector-effect", "non-scaling-stroke");
        brushGroup
            .selectAll(".handle")
            .attr("fill", handleFill)
            .attr("stroke", "none")
            .attr("vector-effect", "non-scaling-stroke");

        return () => {
            svg.select(".brush").remove();
            brushRef.current = null;
        };
    }, [
        handleFill,
        handleSize,
        height,
        minMax,
        ref,
        selectionFill,
        selectionStroke,
        setValue,
        updateOn,
        width,
        x,
        xScaleType,
        y,
    ]);

    const brushValue = brushConfig?.value;

    const setBrushValue = useCallback(
        (value: Range | null | undefined) => {
            if (!brushRef.current || !ref.current || !minMax) return;
            const svg = d3.select(ref.current);
            const xScale = createBrushScale(xScaleType, minMax, [x, x + width]);

            if (!value) {
                // @ts-ignore d3 brush typings are not worth fighting here
                svg.select(".brush").call(brushRef.current.move, null);
                return;
            }
            const [start, end] = value;
            const x0 = xScale(start);
            const x1 = xScale(end);
            // @ts-ignore d3 brush typings are not worth fighting here
            svg.select(".brush").call(brushRef.current.move, [x0, x1]);
        },
        [minMax, ref, width, x, xScaleType],
    );

    useEffect(() => {
        if (isUserBrushingRef.current) return;
        setBrushValue(brushValue);
    }, [brushValue, setBrushValue]);
}
