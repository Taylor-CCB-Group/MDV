import {
    LEGEND_CATEGORICAL_LABEL_START_X,
    LEGEND_CATEGORICAL_ROW_HEIGHT,
} from "@/react/legend/shared/legendConstants";
import type { LegendCategoricalSvgProps } from "@/react/legend/shared/legendTypes";
import {
    LEGEND_CATEGORICAL_LABEL_MAX_WIDTH,
    formatLegendLabel,
    legendCategoricalBodyHeight,
    legendCategoricalContainerHeight,
    legendCategoricalRowY,
} from "@/react/legend/shared/legendUtils";

export {
    getCategoricalLabelMaxWidth,
    LEGEND_CATEGORICAL_LABEL_MAX_WIDTH,
    legendCategoricalContainerHeight,
} from "@/react/legend/shared/legendUtils";
export { LEGEND_CATEGORICAL_WIDTH } from "@/react/legend/shared/legendConstants";
export type {
    LegendCategoricalSvgItem,
    LegendCategoricalSvgProps,
} from "@/react/legend/shared/legendTypes";

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
}: LegendCategoricalSvgProps) {
    const height = legendCategoricalBodyHeight(items.length);

    return (
        <svg
            height={height}
            width="100%"
            className="relative overflow-hidden"
        >
            <g>
                {items.map((item, i) => {
                    const y = legendCategoricalRowY(i);
                    const formatted = formatLegendLabel(
                        item.label,
                        labelMaxWidth,
                    );
                    const isHovered =
                        hoveredKey !== null && hoveredKey === item.key;
                    const isActive =
                        activeKey !== null && activeKey === item.key;
                    const focusedKey = hoveredKey ?? activeKey;
                    const dimmed =
                        focusedKey !== null &&
                        focusedKey !== item.key &&
                        !isActive;
                    const isInteractive = Boolean(onItemHover || onItemClick);
                    return (
                        <g
                            key={item.key}
                            onMouseEnter={
                                onItemHover
                                    ? () => onItemHover(item.key)
                                    : undefined
                            }
                            onMouseLeave={
                                onItemHover
                                    ? () => onItemHover(null)
                                    : undefined
                            }
                            onClick={
                                onItemClick
                                    ? () => onItemClick(item.key)
                                    : undefined
                            }
                            className={isInteractive ? "cursor-pointer" : undefined}
                            style={{
                                opacity: dimmed ? 0.4 : 1,
                            }}
                        >
                            <rect
                                y={Math.max(0, y - 1)}
                                x={0}
                                height={LEGEND_CATEGORICAL_ROW_HEIGHT + 2}
                                width="100%"
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
                                stroke={
                                    isHovered || isActive
                                        ? "currentColor"
                                        : "none"
                                }
                                strokeWidth={isHovered ? 2 : isActive ? 1.5 : 0}
                            />
                            {formatted.truncated ? (
                                <title>{formatted.full}</title>
                            ) : null}
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
