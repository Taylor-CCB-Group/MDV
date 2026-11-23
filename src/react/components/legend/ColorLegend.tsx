import { useCallback, useLayoutEffect, useRef, useState } from "react";
import LegendCategoricalSvg, {
    getCategoricalLabelMaxWidth,
    LEGEND_CATEGORICAL_WIDTH,
    legendCategoricalContainerHeight,
} from "@/react/components/legend/LegendCategoricalSvg";
import LegendContinuousSvg, { DEFAULT_CONTINUOUS_LEGEND_WIDTH } from "@/react/components/legend/LegendContinuousSvg";
import { type LegendFrameSize, useLegendFrame, useLegendPreferredSize } from "@/react/components/legend/LegendFrame";
import type { ColorLegendSpec } from "@/react/legend/color_legend/types";
import {
    LEGEND_BAND_CATEGORICAL_CHROME,
    LEGEND_BAND_COLUMN_WIDTH,
    LEGEND_BAND_MAX_HEIGHT,
} from "@/react/legend/shared/legendConstants";
import { getContinuousLegendContainerHeight, legendCategoricalBodyHeight } from "@/react/legend/shared/legendUtils";

export type ColorLegendProps = {
    spec: ColorLegendSpec;
    activeCategoricalValue?: string | null;
    onCategoricalItemClick?: (value: string) => void;
    activeContinuousRange?: [number, number] | null;
    onContinuousRangeChange?: (range: [number, number] | null) => void;
};

function ColorLegendCategorical({
    label,
    items,
    activeValue = null,
    onItemClick,
    columns = 1,
}: {
    label: string;
    items: { color: string; name: string; value: string }[];
    activeValue?: string | null;
    onItemClick?: (value: string) => void;
    columns?: number;
}) {
    const [hoveredValue, setHoveredValue] = useState<string | null>(null);
    const titleRef = useRef<HTMLDivElement>(null);
    const bodyRef = useRef<HTMLDivElement>(null);
    const [showTitleTooltip, setShowTitleTooltip] = useState(false);
    const [bodyWidth, setBodyWidth] = useState(LEGEND_CATEGORICAL_WIDTH);
    const columnWidth = bodyWidth / columns;
    const labelMaxWidth = getCategoricalLabelMaxWidth(columnWidth);
    const svgItems = items.map((item) => ({
        key: item.value,
        color: item.color,
        label: item.name,
    }));

    const updateTitleTooltip = useCallback(() => {
        const el = titleRef.current;
        if (!el) {
            return;
        }
        setShowTitleTooltip(el.scrollWidth > el.clientWidth);
    }, []);

    useLayoutEffect(() => {
        updateTitleTooltip();
    });

    useLayoutEffect(() => {
        const body = bodyRef.current;
        if (!body) {
            return;
        }
        const updateLayout = () => {
            setBodyWidth(body.clientWidth || LEGEND_CATEGORICAL_WIDTH);
            updateTitleTooltip();
        };
        updateLayout();
        if (typeof ResizeObserver === "undefined") {
            return;
        }
        const observer = new ResizeObserver(updateLayout);
        observer.observe(body);
        if (titleRef.current) {
            observer.observe(titleRef.current);
        }
        return () => observer.disconnect();
    }, [updateTitleTooltip]);

    return (
        <div className="h-full w-full overflow-hidden">
            <div
                ref={titleRef}
                className="legend-title legend-drag-handle overflow-hidden text-ellipsis font-medium"
                title={showTitleTooltip ? label : undefined}
                style={{
                    height: "20px",
                    whiteSpace: "nowrap",
                }}
            >
                {label}
            </div>
            <div
                ref={bodyRef}
                className="legend-body overflow-y-auto overflow-x-hidden w-full"
                style={{
                    height: "calc(100% - 20px)",
                    marginTop: 4,
                }}
            >
                <LegendCategoricalSvg
                    items={svgItems}
                    hoveredKey={hoveredValue}
                    activeKey={activeValue}
                    onItemHover={setHoveredValue}
                    onItemClick={onItemClick}
                    labelMaxWidth={labelMaxWidth}
                    columns={columns}
                    columnWidth={columnWidth}
                />
            </div>
        </div>
    );
}

/** How big the legend wants to be: its usual size floating, or a height that fits the strip. */
function getPreferredSize(spec: ColorLegendSpec, inBand: boolean, columns: number): LegendFrameSize {
    if (spec.kind === "categorical") {
        if (inBand) {
            // the strip sets the width; the height fits the rows the columns leave, up to the cap
            const rows = Math.ceil(spec.items.length / columns);
            return {
                height: Math.min(
                    LEGEND_BAND_MAX_HEIGHT,
                    legendCategoricalBodyHeight(rows) + LEGEND_BAND_CATEGORICAL_CHROME,
                ),
            };
        }
        return { width: LEGEND_CATEGORICAL_WIDTH, height: legendCategoricalContainerHeight(spec.items.length) };
    }
    if (inBand) {
        // the strip's frame adds 2px padding top and bottom
        const compact = getContinuousLegendContainerHeight(spec.range, {
            hasLabel: Boolean(spec.label),
            compact: true,
        });
        return { height: compact + 4 };
    }
    const width = spec.width ?? DEFAULT_CONTINUOUS_LEGEND_WIDTH;
    const height =
        spec.height ??
        getContinuousLegendContainerHeight(spec.range, {
            width,
            hasLabel: Boolean(spec.label),
            isDate: spec.isDate,
        });
    return { width, height };
}

/**
 * Inner content for the chart color legend (categorical list or continuous color bar). In the
 * strip above the plot it fills the width it's given, with a categorical list in columns.
 */
export default function ColorLegend({
    spec,
    activeCategoricalValue = null,
    onCategoricalItemClick,
    activeContinuousRange = null,
    onContinuousRangeChange,
}: ColorLegendProps) {
    const { inBand, width: bandWidth } = useLegendFrame();
    const itemCount = spec.kind === "categorical" ? spec.items.length : 0;
    const columns = inBand ? Math.max(1, Math.min(itemCount, Math.floor(bandWidth / LEGEND_BAND_COLUMN_WIDTH))) : 1;
    useLegendPreferredSize(getPreferredSize(spec, inBand, columns));

    if (spec.kind === "categorical") {
        return (
            <div className="h-full w-full">
                <ColorLegendCategorical
                    label={spec.label}
                    items={spec.items}
                    activeValue={activeCategoricalValue}
                    onItemClick={onCategoricalItemClick}
                    columns={columns}
                />
            </div>
        );
    }

    return (
        <div className="h-full w-full overflow-hidden">
            <LegendContinuousSvg
                label={spec.label}
                colors={spec.colors}
                range={spec.range}
                width={spec.width}
                height={spec.height}
                activeRange={activeContinuousRange}
                onRangeChange={onContinuousRangeChange}
                isDate={spec.isDate}
                compact={inBand}
            />
        </div>
    );
}
