import { observer } from "mobx-react-lite";
import { ScatterPlotConfig2D, ScatterPlotConfig3D } from "../scatter_state";
import type { DataColumn } from "@/charts/charts";
import { useChartSize, useParamColumns } from "../hooks";
import { useMemo, type PropsWithChildren } from "react";
import { Axis, Scale } from "@visx/visx";

type AxisComponentProps = {
    config: ScatterPlotConfig2D | ScatterPlotConfig3D
    unproject: (xy: [number, number]) => number[];
} & PropsWithChildren;

export default observer(function AxisComponent({ config, unproject, children }: AxisComponentProps) {
    const [cx, cy] = useParamColumns() as DataColumn<"double">[];
    const { dimension, viewState } = config;
    const is2d = dimension === "2d";
    const margin = useMemo(() => (
        //todo better Axis encapsulation in another component
        is2d ? {
            top: 10,
            right: 10,
            bottom: config.axis.x.size + 20,
            left: config.axis.y.size + 20,
        } : {
            top: 0,
            right: 0,
            bottom: 0,
            left: 0,
        }
    ), [is2d, is2d && config.axis.x.size, is2d && config.axis.y.size]);
    const [width, height] = useChartSize();
    const chartWidth = width - margin.left - margin.right;
    //there could be a potential off-by-one/two error somewhere down the line
    //if we don't fully understand reasons for `- 2` here.
    //prevents overlapping with x-axis.
    const chartHeight = height - margin.top - margin.bottom - 2;

    // axes need to respond to the viewState... (make sure there isn't a regression here when refactoring etc).
    const ranges = useMemo(() => {
        viewState;
        // first time around, we get an exception because scatterplotLayer hasn't been rendered yet
        try {
            const p = unproject([0, 0]);
            const p2 = unproject([chartWidth, chartHeight]);
            const domainX = [p[0], p2[0]];
            const domainY = [p2[1], p[1]];
            return { domainX, domainY };
        } catch (e) {
            return { domainX: cx.minMax, domainY: cy.minMax };
        }
    }, [cx.minMax, cy.minMax, viewState, chartWidth, chartHeight, unproject]);
    // * as of now, we only use these scales for the axes,
    // but we should consider how they might relate to data transformation
    const scaleX = useMemo(() => Scale.scaleLinear({
        domain: ranges.domainX, // e.g. [min, max]
        range: [margin.left, chartWidth + margin.left],
    }), [chartWidth, ranges]);
    const scaleY = useMemo(() => Scale.scaleLinear({
        domain: ranges.domainY, // e.g. [min, max]
        range: [chartHeight + margin.top, margin.top],
    }), [chartHeight, ranges]);

    const deckStyle = useMemo(() => ({
        position: "absolute",
        top: margin.top,
        left: margin.left,
        width: chartWidth,
        height: chartHeight,
    } as const), [chartWidth, chartHeight]);
    return (
        <>
            <div style={deckStyle}>{children}</div>
            {is2d && <svg width={width} height={height}>
                <Axis.AxisBottom
                    top={chartHeight + margin.top}
                    scale={scaleX}
                    stroke={"var(--text_color)"}
                    tickStroke={"var(--text_color)"}
                    tickLabelProps={() => ({
                        fill: "var(--text_color)",
                        fontSize: config.axis.x.tickfont,
                        textAnchor: "middle",
                    })}
                    labelProps={{
                        fill: "var(--text_color)",
                        fontSize: config.axis.x.tickfont,
                        textAnchor: "middle",
                    }}
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
                    })}
                    labelProps={{
                        fill: "var(--text_color)",
                        fontSize: config.axis.y.tickfont,
                        // dx: config.axis.x.rotate_labels ? 
                    }}
                    label={cy.name}
                />
            </svg>}        
        </>
    )
});