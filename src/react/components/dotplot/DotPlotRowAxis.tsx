import { AxisLeft } from "@visx/axis";
import type { ScaleBand } from "d3-scale";

type DotPlotRowAxisProps = {
    /** Row category codes, top to bottom. */
    rows: number[];
    /** Label for each row category code. */
    labels: readonly string[];
    scale: ScaleBand<number>;
    left: number;
    top: number;
    /** Label font size; labels shrink to the row height if rows are thinner. */
    tickfont: number;
};

/** The row labels down the left of a dot plot: one per row, never thinned out. */
export default function DotPlotRowAxis({ rows, labels, scale, left, top, tickfont }: DotPlotRowAxisProps) {
    return (
        <AxisLeft
            left={left}
            top={top}
            scale={scale}
            // Without tickValues visx aims for about 10 ticks and drops labels from
            // band scales with more rows than that.
            tickValues={rows}
            tickFormat={(row) => labels[row]}
            stroke="var(--text_color)"
            tickStroke="var(--text_color)"
            tickLabelProps={{
                fill: "var(--text_color)",
                // shrink labels to the row height so they don't overlap on short charts
                fontSize: Math.min(tickfont, scale.bandwidth()),
                textAnchor: "end",
            }}
        />
    );
}
