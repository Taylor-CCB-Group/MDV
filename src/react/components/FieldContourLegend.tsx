import { useEffect, useId, useRef, useState } from "react";
import type { FieldName } from "@/charts/charts";
import type { HoveredContour } from "../contour_state";

export interface FieldLegendItem {
    name: string;
    field: FieldName;
    color: [number, number, number];
}

export interface CategoryLegendItem {
    /** id of the categorical contour layer this item describes */
    id: string;
    name: string;
    /** a single colour, or the colour range shared by a contour over several categories */
    colors: readonly (readonly number[])[];
}

export interface FieldContourLegendProps {
    fields: FieldLegendItem[];
    categories?: CategoryLegendItem[];
    /** heading for the categorical contours, typically the contour parameter's name */
    categoryLabel?: string;
    label?: string;
    position?: { x: number; y: number };
    onContourHover?: (contour: HoveredContour | null) => void;
}

type LegendRow =
    | { type: "heading"; key: string; text: string }
    | { type: "category"; key: string; item: CategoryLegendItem }
    | { type: "field"; key: string; item: FieldLegendItem };

function isSameContour(a: HoveredContour | null, b: HoveredContour) {
    if (a?.type === "field" && b.type === "field") return a.field === b.field;
    if (a?.type === "category" && b.type === "category") return a.id === b.id;
    return false;
}

// Convert RGB array to CSS color string
const rgbToCss = (rgb: readonly number[]): string => {
    return `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
};

/**
 * React component for displaying a legend of density contours with their colors:
 * categorical contours (grouped under the contour parameter's name) and density fields.
 * Hovering an item highlights its contour via `onContourHover`.
 * The legend is fixed in position (typically bottom-left).
 *
 * Note: SVG is used for the legend content to ensure compatibility with figure export
 * functionality, allowing legends to be included in exported visualizations.
 */
export default function FieldContourLegend({
    fields,
    categories = [],
    categoryLabel = "Categories",
    label = categories.length > 0 ? "Density" : "Density Fields",
    position,
    onContourHover,
}: FieldContourLegendProps) {
    const containerRef = useRef<HTMLDivElement>(null);
    const bodyRef = useRef<HTMLDivElement>(null);
    const [hoveredContour, setHoveredContour] = useState<HoveredContour | null>(null);
    // gradient ids must be unique within the document; useId's colons are awkward in url(#...)
    const gradientIdPrefix = `contour-legend-${useId().replace(/:/g, "")}`;

    useEffect(() => {
        const container = containerRef.current;
        if (!container) return;

        // Set fixed position if provided
        if (position) {
            container.style.left = `${position.x}px`;
            container.style.bottom = `${position.y}px`;
        }
    }, [position]);

    if (fields.length === 0 && categories.length === 0) {
        return null;
    }

    const rows: LegendRow[] = [];
    if (categories.length > 0) {
        rows.push({ type: "heading", key: "heading-categories", text: categoryLabel });
        for (const item of categories) rows.push({ type: "category", key: item.id, item });
        // only distinguish fields with their own heading when they share the legend with categories
        if (fields.length > 0) rows.push({ type: "heading", key: "heading-fields", text: "Fields" });
    }
    for (const item of fields) rows.push({ type: "field", key: item.field, item });

    const len = rows.length;
    const h_fac = 12; // Height factor for each item
    const height = len * h_fac + (len + 2) * 2;
    const maxHeight = 215;
    const containerHeight = height > maxHeight ? maxHeight + 35 : height + 35;

    const handleMouseEnter = (contour: HoveredContour) => {
        setHoveredContour(contour);
        onContourHover?.(contour);
    };

    const handleMouseLeave = () => {
        setHoveredContour(null);
        onContourHover?.(null);
    };

    // props for the row and swatch of an item whose contour can be hovered
    const getHoverProps = (contour: HoveredContour) => {
        const isHovered = isSameContour(hoveredContour, contour);
        return {
            groupProps: {
                onMouseEnter: () => handleMouseEnter(contour),
                onMouseLeave: handleMouseLeave,
                className: "cursor-pointer",
                style: {
                    opacity: isHovered ? 1 : hoveredContour ? 0.4 : 1,
                },
            },
            swatchProps: {
                stroke: isHovered ? "currentColor" : "none",
                strokeWidth: isHovered ? 2 : 0,
            },
        };
    };

    const renderRow = (row: LegendRow, y: number) => {
        if (row.type === "heading") {
            return (
                <text
                    key={row.key}
                    y={y + 6 + 3}
                    x={2}
                    className="text-xs fill-current font-bold"
                >
                    {row.text}
                </text>
            );
        }
        if (row.type === "category") {
            const { item } = row;
            const gradientId = `${gradientIdPrefix}-${item.id}`;
            const isGradient = item.colors.length > 1;
            const { groupProps, swatchProps } = getHoverProps({ type: "category", id: item.id });
            return (
                <g key={row.key} {...groupProps}>
                    <title>{item.name}</title>
                    {isGradient && (
                        <defs>
                            <linearGradient id={gradientId}>
                                {item.colors.map((color, i) => (
                                    <stop
                                        // biome-ignore lint/suspicious/noArrayIndexKey: stops are positional
                                        key={i}
                                        offset={i / (item.colors.length - 1)}
                                        stopColor={rgbToCss(color)}
                                    />
                                ))}
                            </linearGradient>
                        </defs>
                    )}
                    <rect
                        y={y}
                        x={2}
                        height="10"
                        width="10"
                        fill={isGradient ? `url(#${gradientId})` : rgbToCss(item.colors[0])}
                        {...swatchProps}
                    />
                    <text
                        y={y + 6 + 3}
                        x={14}
                        className="text-xs fill-current"
                    >
                        {item.name}
                    </text>
                </g>
            );
        }
        const field = row.item;
        const color = rgbToCss(field.color);
        const { groupProps, swatchProps } = getHoverProps({ type: "field", field: field.field });
        return (
            <g key={row.key} {...groupProps}>
                <rect
                    y={y}
                    x={2}
                    height="10"
                    width="10"
                    fill={color}
                    {...swatchProps}
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
                className="legend-title h-5 whitespace-nowrap px-1 text-xs font-bold text-current"
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
                        {rows.map((row, i) => renderRow(row, (i + 1) * 2 + i * h_fac))}
                    </g>
                </svg>
            </div>
        </div>
    );
}
