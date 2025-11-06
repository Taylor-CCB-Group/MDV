import { observer } from "mobx-react-lite";
import type { AxisConfig, ScatterPlotConfig2D, ScatterPlotConfig3D } from "../scatter_state";
import type { DataColumn } from "@/charts/charts";
import { useChartSize, useParamColumns } from "../hooks";
import { useEffect, useLayoutEffect, useMemo, useState, type PropsWithChildren } from "react";
// import * as Axis from "@visx/axis";
// import * as Scale from "@visx/scale";

//! This component is broken as the visx dependency is uninstalled due to incompatibility with React 19 version
type AxisComponentProps = {
    config: ScatterPlotConfig2D | ScatterPlotConfig3D
    unproject: (xy: [number, number]) => number[];
} & PropsWithChildren;

//! todo - get label options working... for now, we ignore rotate_labels
// not generally relevant to scatterplot, but will be for other charts
function getLabelProps(axisConfig?: AxisConfig) {
    if (!axisConfig) return {};
    const { tickfont, rotate_labels } = axisConfig;
    return () => ({
        fill: "var(--text_color)",
        fontSize: tickfont,
        textAnchor: rotate_labels ? "end" : "middle",
        dx: rotate_labels ? "-.8em" : undefined,
        dy: rotate_labels ? "-.10em" : ".71em",
        transform: rotate_labels ? "rotate(-45)" : undefined,
    });
}

// todo export this & avoid duplication in DeckScatterComponent, need a little bit of a refactor for right props
function useSynchronizedScales({ config, unproject }: AxisComponentProps) {
    const [cx, cy] = useParamColumns() as DataColumn<"double">[];
    const { dimension, viewState } = config;
    const is2d = dimension === "2d";
    const xSize = is2d ? config.axis.x.size : 0;
    const ySize = is2d ? config.axis.y.size : 0;

    const margin = useMemo(() => (
        is2d ? {
            top: 10,
            right: 10,
            bottom: xSize,
            left: ySize,
        } : {
            top: 0,
            right: 0,
            bottom: 0,
            left: 0,
        }
    ), [is2d, xSize, ySize]);
    const [width, height] = useChartSize();
    const chartWidth = width - margin.left - margin.right;
    //there could be a potential off-by-one/two error somewhere down the line
    //if we don't fully understand reasons for `- 2` here.
    //prevents overlapping with x-axis.
    const chartHeight = height - margin.top - margin.bottom - 2;

    const [ranges, setRanges] = useState({ domainX: cx.minMax, domainY: cy.minMax });

    // Run synchronously after layout/before paint
    useLayoutEffect(() => {
        //! not reacting to changes from useZoomOnFilter()...
        viewState;
        // first time around, we get an exception because scatterplotLayer hasn't been rendered yet
        try {
            const p = unproject([0, 0]);
            const p2 = unproject([chartWidth, chartHeight]);
            setRanges({
                domainX: [p[0], p2[0]],
                domainY: [p2[1], p[1]],
            });
        } catch (e) {
            // Fallback to data ranges
            // console.warn("AxisComponent: unproject failed", e);
            setRanges({ domainX: cx.minMax, domainY: cy.minMax });
        }
    }, [viewState, chartWidth, chartHeight, unproject, cx.minMax, cy.minMax]);

    // const scaleX = useMemo(() => Scale.scaleLinear({
    //     domain: ranges.domainX.sort(), // e.g. [min, max]
    //     range: [margin.left, chartWidth + margin.left],
    // }), [chartWidth, ranges, margin.left]);
    // const scaleY = useMemo(() => Scale.scaleLinear({
    //     domain: ranges.domainY.sort(), // e.g. [min, max]
    //     range: [chartHeight + margin.top, margin.top],
    // }), [chartHeight, ranges, margin.top]);

    return { 
        ranges, margin, chartWidth, chartHeight, 
        // scaleX, scaleY 
    };
}

export default observer(function AxisComponent({ config, unproject, children }: AxisComponentProps) {
    const [cx, cy] = useParamColumns() as DataColumn<"double">[];
    const { dimension } = config;
    const is2d = dimension === "2d";
    const [width, height] = useChartSize();
    const { 
        // scaleX, scaleY, 
        margin, chartWidth, chartHeight } = useSynchronizedScales({ config, unproject });
    
    const deckStyle = useMemo(() => ({
        position: "absolute",
        top: margin.top,
        left: margin.left,
        width: chartWidth,
        height: chartHeight,
    } as const), [chartWidth, chartHeight, margin.top, margin.left]);
    // useEffect(() => {
    //     if (is2d && (config.axis.x.rotate_labels || config.axis.y.rotate_labels)) {
    //         console.warn("Axis rotation not implemented for react charts");
    //     }
    // }, [is2d]);
    useEffect(() => {
        console.warn("This is a broken chart, visx dependency has been removed");
    }, []);
    return (
        <>
            <div style={deckStyle}>{children}</div>
            {is2d && <svg width={width} height={height}>
                {/* <Axis.AxisBottom
                    top={chartHeight + margin.top}
                    scale={scaleX}
                    stroke={"var(--text_color)"}
                    tickStroke={"var(--text_color)"}
                    tickLabelProps={() => ({
                        fill: "var(--text_color)",
                        fontSize: config.axis.x.tickfont,
                        textAnchor: "middle",
                        // for some reason, tickfont triggers a re-render but rotate_labels doesn't
                        angle: config.axis.x.rotate_labels ? -45 : 0,
                    })}
                    labelProps={{
                        fill: "var(--text_color)",
                        fontSize: config.axis.x.tickfont,
                        textAnchor: "middle",
                    }}
                    labelOffset={0}
                    label={cx.name}
                />
                <Axis.AxisLeft
                    left={margin.left}
                    scale={scaleY}
                    stroke={"var(--text_color)"}
                    tickStroke={"var(--text_color)"}
                    tickLabelProps={() => ({
                        fill: "var(--text_color)",
                        fontSize: config.axis.y.tickfont,
                        textAnchor: "end",
                        angle: config.axis.y.rotate_labels ? -45 : 0,
                    })}
                    labelProps={{
                        fill: "var(--text_color)",
                        fontSize: config.axis.y.tickfont,
                    }}
                    labelOffset={20}
                    label={cy.name}
                /> */}
            </svg>}        
        </>
    )
});