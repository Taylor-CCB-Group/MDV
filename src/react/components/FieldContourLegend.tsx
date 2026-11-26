import { useId, useMemo, useState } from "react";
import type { FieldName } from "@/charts/charts";
import { type LegendSummary, useLegendFrame } from "@/react/components/legend/LegendFrame";
import type { LegendWrapperComponentProps } from "@/react/components/legend/LegendWrapper";
import { type ChartLegendConfig, useChartLegend } from "@/react/components/legend/useChartLegend";
import { useChart } from "@/react/context";
import { LEGEND_BAND_MAX_HEIGHT } from "@/react/legend/shared/legendConstants";
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
 * Charts show it with `useContourLegend`, in the shared legend wrapper (draggable, can be minimised).
 *
 * Note: SVG is used for the legend content to ensure compatibility with figure export
 * functionality, allowing legends to be included in exported visualizations.
 */
export default function FieldContourLegend({
    fields,
    categories = [],
    categoryLabel = "Categories",
    label = categories.length > 0 ? "Density" : "Density Fields",
    onContourHover,
}: FieldContourLegendProps) {
    const [hoveredContour, setHoveredContour] = useState<HoveredContour | null>(null);
    // gradient ids must be unique within the document; useId's colons are awkward in url(#...)
    const gradientIdPrefix = `contour-legend-${useId().replace(/:/g, "")}`;

    const { inBand } = useLegendFrame();

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
                <text key={row.key} y={y + 6 + 3} x={2} className="text-xs fill-current font-bold">
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
                    <text y={y + 6 + 3} x={14} className="text-xs fill-current">
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
                <rect y={y} x={2} height="10" width="10" fill={color} {...swatchProps} />
                <text y={y + 6 + 3} x={14} className="text-xs fill-current">
                    {field.name === "" ? "none" : field.name}
                </text>
            </g>
        );
    };

    return (
        <div
            style={{
                width: inBand ? "100%" : "140px",
                // the strip's wrapper adds 2px padding top and bottom
                height: `${inBand ? Math.min(containerHeight, LEGEND_BAND_MAX_HEIGHT - 4) : containerHeight}px`,
            }}
        >
            <div className="legend-title legend-drag-handle h-5 whitespace-nowrap pl-1 pr-12 text-xs font-bold text-current">
                {label}
            </div>
            <div
                className="legend-body overflow-y-auto overflow-x-hidden w-full"
                style={{
                    height: "calc(100% - 25px)",
                }}
            >
                <svg height={height} width={inBand ? "100%" : 180} className="relative">
                    <g>{rows.map((row, i) => renderRow(row, (i + 1) * 2 + i * h_fac))}</g>
                </svg>
            </div>
        </div>
    );
}

export type ContourLegendSpec = FieldContourLegendProps;

/** The contour legend as the legend wrapper draws it. */
function ContourLegendInWrapper({ spec }: LegendWrapperComponentProps<ContourLegendSpec>) {
    return <FieldContourLegend {...spec} />;
}

/** The minimised contour legend: its title and a few of the contours' colours. */
function contourLegendSummary(spec: ContourLegendSpec): LegendSummary {
    const colors = [
        ...(spec.categories ?? []).map((item) => rgbToCss(item.colors[0])),
        ...spec.fields.map((item) => rgbToCss(item.color)),
    ];
    const shown = colors.slice(0, 4);
    return {
        label: spec.label ?? ((spec.categories?.length ?? 0) > 0 ? "Density" : "Density Fields"),
        swatches: shown,
        more: colors.length - shown.length,
    };
}

/**
 * Shows a chart's contour legend in the shared legend wrapper: draggable, can be minimised, and
 * goes in the strip above the plot with the chart's other legends. Starts at the bottom left.
 */
export function useContourLegend(
    legend: ChartLegendConfig,
    { fields, categories = [], categoryLabel, onContourHover }: ContourLegendSpec,
) {
    const chart = useChart();
    const spec = useMemo<ContourLegendSpec | null>(
        () =>
            fields.length > 0 || categories.length > 0 ? { fields, categories, categoryLabel, onContourHover } : null,
        [fields, categories, categoryLabel, onContourHover],
    );
    const summary = useMemo<LegendSummary>(() => (spec ? contourLegendSummary(spec) : { label: "Density" }), [spec]);
    useChartLegend({
        chart,
        legend,
        spec,
        component: ContourLegendInWrapper,
        summary,
        defaultBottomLeft: [10, 10],
    });
}
