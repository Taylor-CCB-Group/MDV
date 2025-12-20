import { useEffect, useRef, useState } from "react";
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
    onPositionChange?: (position: { x: number; y: number }) => void;
    onFieldHover?: (fieldId: FieldName | null) => void;
}

/**
 * React component for displaying a legend of field contours with their colors.
 * The legend is draggable and resizable, similar to the existing getColorLegend utility.
 * 
 * Note: SVG is used for the legend content to ensure compatibility with figure export
 * functionality, allowing legends to be included in exported visualizations.
 */
export default function FieldContourLegend({
    fields,
    label = "Density Fields",
    position,
    onPositionChange,
    onFieldHover,
}: FieldContourLegendProps) {
    const containerRef = useRef<HTMLDivElement>(null);
    const bodyRef = useRef<HTMLDivElement>(null);
    const [hoveredField, setHoveredField] = useState<FieldName | null>(null);

    useEffect(() => {
        const container = containerRef.current;
        if (!container) return;

        // Set initial position if provided
        if (position) {
            container.style.left = `${position.x}px`;
            container.style.top = `${position.y}px`;
        }

        // Apply draggable and resizable functionality
        // Using the same utilities as the existing legend system
        // Header (legend-title) is draggable for better UX
        const dragConfig: any = { handle: ".legend-title" };
        
        // Track position changes when dragging ends
        if (onPositionChange) {
            dragConfig.ondragend = () => {
                // Get position from the container element
                const currentContainer = containerRef.current;
                if (currentContainer) {
                    onPositionChange({
                        x: currentContainer.offsetLeft,
                        y: currentContainer.offsetTop,
                    });
                }
            };
        }
        
        makeDraggable(container, dragConfig);
        makeResizable(container);
        
        // Cleanup function to remove event listeners if needed
        return () => {
            // makeDraggable doesn't provide a cleanup function, but the element
            // will be cleaned up when the component unmounts
        };
    }, [position, onPositionChange]);

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

    const handleMouseEnter = (fieldId: FieldName) => {
        setHoveredField(fieldId);
        onFieldHover?.(fieldId);
    };

    const handleMouseLeave = () => {
        setHoveredField(null);
        onFieldHover?.(null);
    };

    return (
        <div
            ref={containerRef}
            className="legend-container absolute border-[0.5px] border-current z-[2]"
            style={{
                width: "120px",
                height: `${containerHeight}px`,
            }}
        >
            <div
                className="legend-title h-5 whitespace-nowrap px-1 text-xs font-bold text-current cursor-move"
            >
                {label}
            </div>
            <div
                ref={bodyRef}
                className="legend-body overflow-y-auto overflow-x-hidden w-full"
                style={{
                    height: "calc(100% - 25px)",
                }}
            >
                <svg
                    height={height}
                    width={180}
                    className="relative"
                >
                    <g>
                        {fields.map((field, i) => {
                            const y = (i + 1) * 2 + i * h_fac;
                            const color = rgbToCss(field.color);
                            const isHovered = hoveredField === field.field;
                            return (
                                <g 
                                    key={field.field}
                                    onMouseEnter={() => handleMouseEnter(field.field)}
                                    onMouseLeave={handleMouseLeave}
                                    className="cursor-pointer"
                                    style={{
                                        opacity: isHovered ? 1 : hoveredField ? 0.4 : 1,
                                    }}
                                >
                                    <rect
                                        y={y}
                                        x={2}
                                        height="10"
                                        width="10"
                                        fill={color}
                                        stroke={isHovered ? "currentColor" : "none"}
                                        strokeWidth={isHovered ? 2 : 0}
                                    />
                                    <text
                                        y={y + 6 + 3}
                                        x={14}
                                        className="text-xs fill-current"
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

