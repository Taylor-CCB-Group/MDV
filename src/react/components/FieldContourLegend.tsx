import { useEffect, useRef } from "react";
import type { FieldName } from "@/charts/charts";
import { makeDraggable, makeResizable } from "@/utilities/Elements";

export interface FieldLegendItem {
    name: string;
    field: FieldName;
    color: [number, number, number];
}

export interface FieldContourLegendProps {
    fields: FieldLegendItem[];
    label?: string;
    position?: { x: number; y: number };
}

/**
 * React component for displaying a legend of field contours with their colors.
 * The legend is draggable and resizable, similar to the existing getColorLegend utility.
 */
export default function FieldContourLegend({
    fields,
    label = "Density Fields",
    position,
}: FieldContourLegendProps) {
    const containerRef = useRef<HTMLDivElement>(null);
    const bodyRef = useRef<HTMLDivElement>(null);

    useEffect(() => {
        const container = containerRef.current;
        if (!container) return;

        // Apply draggable and resizable functionality
        // Using the same utilities as the existing legend system
        makeDraggable(container, { handle: ".legend-body" });
        makeResizable(container);

        // Set initial position if provided
        if (position) {
            container.style.left = `${position.x}px`;
            container.style.top = `${position.y}px`;
        }
    }, [position]);

    if (fields.length === 0) {
        return null;
    }

    const len = fields.length;
    const h_fac = 12; // Height factor for each item
    const height = len * h_fac + (len + 2) * 2;
    const maxHeight = 215;
    const containerHeight = height > maxHeight ? maxHeight + 35 : height + 35;

    // Convert RGB array to CSS color string
    const rgbToCss = (rgb: [number, number, number]): string => {
        return `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
    };

    return (
        <div
            ref={containerRef}
            className="legend-container"
            style={{
                width: "120px",
                height: `${containerHeight}px`,
                position: "absolute",
                border: "0.5px solid currentcolor",
                zIndex: 2,
            }}
        >
            <div
                className="legend-title"
                style={{
                    height: "20px",
                    whiteSpace: "nowrap",
                    padding: "2px 4px",
                    fontSize: "12px",
                    fontWeight: "bold",
                    color: "currentcolor",
                }}
            >
                {label}
            </div>
            <div
                ref={bodyRef}
                className="legend-body"
                style={{
                    overflowY: "auto",
                    overflowX: "hidden",
                    height: "calc(100% - 25px)",
                    width: "100%",
                }}
            >
                <svg
                    height={height}
                    width={180}
                    style={{
                        position: "relative",
                    }}
                >
                    <g>
                        {fields.map((field, i) => {
                            const y = (i + 1) * 2 + i * h_fac;
                            const color = rgbToCss(field.color);
                            return (
                                <g key={field.field}>
                                    <rect
                                        y={y}
                                        x={2}
                                        height="10"
                                        width="10"
                                        fill={color}
                                    />
                                    <text
                                        y={y + 6 + 3}
                                        x={14}
                                        alignmentBaseline="central"
                                        style={{
                                            fontSize: "12px",
                                            fill: "currentcolor",
                                        }}
                                    >
                                        {field.name === "" ? "none" : field.name}
                                    </text>
                                </g>
                            );
                        })}
                    </g>
                </svg>
            </div>
        </div>
    );
}

