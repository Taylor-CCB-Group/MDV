import React from "react";
import { ScaleLinear, scaleLinear } from "d3-scale";

// Helper to generate ticks for a linear scale
function getTicks(scale: ScaleLinear<number, number, never>, count = 5) {
  return scale.ticks ? scale.ticks(count) : [];
}

interface SVGChartWrapperProps {
  width?: number;
  data: any; // Assuming data is an array of objects
  height?: number;
  margin?: { top: number; right: number; bottom: number; left: number };
  showXAxis?: boolean;
  showYAxis?: boolean;
  xDomain?: [number, number];
  yDomain?: [number, number];
  children?: React.ReactNode | ((params: {
    xScale: ScaleLinear<number, number, never>;
    yScale: ScaleLinear<number, number, never>;
    innerWidth: number;
    innerHeight: number;
    data:any;
  }) => React.ReactNode);
}

export default function SVGChartWrapper({
  width = 400,
  height = 300,
  data = null,
  margin = { top: 30, right: 30, bottom: 40, left: 40 },
  showXAxis = true,
  showYAxis = true,
  xDomain = [0, 10],
  yDomain = [0, 100],
  children,
}: SVGChartWrapperProps) {
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;

  const xScale = scaleLinear().domain(xDomain).range([0, innerWidth]);
  const yScale = scaleLinear().domain(yDomain).range([innerHeight, 0]);

  const xTicks = getTicks(xScale, 5);
  const yTicks = getTicks(yScale, 5);

  return (
    <svg width={width} height={height} style={{ background: "#fff" }}>
      <g transform={`translate(${margin.left},${margin.top})`}>
        {/* X Axis */}
        {showXAxis && (
          <>
            <line
              x1={0}
              y1={innerHeight}
              x2={innerWidth}
              y2={innerHeight}
              stroke="black"
            />
            {xTicks.map((tick, i) => (
              <g key={i} transform={`translate(${xScale(tick)},${innerHeight})`}>
                <line y2="6" stroke="black" />
                <text y="20" textAnchor="middle" fontSize="10">
                  {tick}
                </text>
              </g>
            ))}
          </>
        )}
        {/* Y Axis */}
        {showYAxis && (
          <>
            <line x1={0} y1={0} x2={0} y2={innerHeight} stroke="black" />
            {yTicks.map((tick, i) => (
              <g key={i} transform={`translate(0,${yScale(tick)})`}>
                <line x2="-6" stroke="black" />
                <text x="-10" y="4" textAnchor="end" fontSize="10">
                  {tick}
                </text>
              </g>
            ))}
          </>
        )}
        {/* Chart content */}
        {typeof children === "function"
          ? children({ xScale, yScale, innerWidth, innerHeight,data })
          : children}
      </g>
    </svg>
  );
}