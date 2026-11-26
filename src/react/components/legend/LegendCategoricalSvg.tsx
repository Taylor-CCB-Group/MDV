import { LEGEND_CATEGORICAL_LABEL_START_X, LEGEND_CATEGORICAL_ROW_HEIGHT } from "@/react/legend/shared/legendConstants";
import type { LegendCategoricalSvgProps } from "@/react/legend/shared/legendTypes";
import {
    formatLegendLabel,
    LEGEND_CATEGORICAL_LABEL_MAX_WIDTH,
    legendCategoricalBodyHeight,
    legendCategoricalRowY,
} from "@/react/legend/shared/legendUtils";

export { LEGEND_CATEGORICAL_WIDTH } from "@/react/legend/shared/legendConstants";
export type {
    LegendCategoricalSvgItem,
    LegendCategoricalSvgProps,
} from "@/react/legend/shared/legendTypes";
export {
    getCategoricalLabelMaxWidth,
    LEGEND_CATEGORICAL_LABEL_MAX_WIDTH,
    legendCategoricalContainerHeight,
} from "@/react/legend/shared/legendUtils";

/**
 * Shared SVG list for categorical legends (color legend and field contour legend).
 */
export default function LegendCategoricalSvg({
    items,
    onItemHover,
    onItemClick,
    hoveredKey = null,
    activeKey = null,
    labelMaxWidth = LEGEND_CATEGORICAL_LABEL_MAX_WIDTH,
    columns = 1,
    columnWidth = 0,
}: LegendCategoricalSvgProps) {
    // filled down each column in turn, so the list still reads top to bottom
    const rows = Math.ceil(items.length / Math.max(1, columns));
    const height = legendCategoricalBodyHeight(rows);

    return (
        <svg height={height} width="100%" className="relative overflow-hidden">
            <g>
                {items.map((item, i) => {
                    const y = legendCategoricalRowY(i % rows);
                    const x = columns > 1 ? Math.floor(i / rows) * columnWidth : 0;
                    const formatted = formatLegendLabel(item.label, labelMaxWidth);
                    const isHovered = hoveredKey !== null && hoveredKey === item.key;
                    const isActive = activeKey !== null && activeKey === item.key;
                    const focusedKey = hoveredKey ?? activeKey;
                    const dimmed = focusedKey !== null && focusedKey !== item.key && !isActive;
                    const isInteractive = Boolean(onItemHover || onItemClick);
                    return (
                        // biome-ignore lint/a11y/noStaticElementInteractions: pointer-only legend rows
                        <g
                            key={item.key}
                            onMouseEnter={onItemHover ? () => onItemHover(item.key) : undefined}
                            onMouseLeave={onItemHover ? () => onItemHover(null) : undefined}
                            onClick={onItemClick ? () => onItemClick(item.key) : undefined}
                            className={isInteractive ? "cursor-pointer" : undefined}
                            transform={x ? `translate(${x}, 0)` : undefined}
                            style={{
                                opacity: dimmed ? 0.4 : 1,
                            }}
                        >
                            <rect
                                y={Math.max(0, y - 1)}
                                x={0}
                                height={LEGEND_CATEGORICAL_ROW_HEIGHT + 2}
                                width={columns > 1 ? columnWidth : "100%"}
                                style={{
                                    fill: "transparent",
                                    pointerEvents: "all",
                                }}
                            />
                            <rect
                                y={y}
                                x={2}
                                height="10"
                                width="10"
                                fill={item.color}
                                stroke={isHovered || isActive ? "currentColor" : "none"}
                                strokeWidth={isHovered ? 2 : isActive ? 1.5 : 0}
                            />
                            {formatted.truncated ? <title>{formatted.full}</title> : null}
                            <text
                                y={y + 9}
                                x={LEGEND_CATEGORICAL_LABEL_START_X}
                                aria-label={formatted.full}
                                className="text-xs fill-current"
                                style={{
                                    fontWeight: isActive ? 500 : undefined,
                                }}
                            >
                                {formatted.display}
                            </text>
                        </g>
                    );
                })}
            </g>
        </svg>
    );
}
