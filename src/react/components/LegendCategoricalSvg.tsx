export type LegendCategoricalSvgItem = {
    key: string;
    color: string;
    label: string;
};

export type LegendCategoricalSvgProps = {
    items: LegendCategoricalSvgItem[];
    /** Called when pointer enters/leaves a row (optional, used by field contour legend). */
    onItemHover?: (key: string | null) => void;
    onItemClick?: (key: string) => void;
    hoveredKey?: string | null;
    activeKey?: string | null;
};

const H_FAC = 12;
const MAX_HEIGHT = 215;
export const LEGEND_CATEGORICAL_WIDTH = 180;

export function legendCategoricalBodyHeight(itemCount: number): number {
    return itemCount * H_FAC + (itemCount + 2) * 2;
}

export function legendCategoricalContainerHeight(itemCount: number): number {
    const bodyHeight = legendCategoricalBodyHeight(itemCount);
    return bodyHeight > MAX_HEIGHT ? MAX_HEIGHT + 35 : bodyHeight + 35;
}

/**
 * Shared SVG list for categorical legends (color legend and field contour legend).
 */
export default function LegendCategoricalSvg({
    items,
    onItemHover,
    onItemClick,
    hoveredKey = null,
    activeKey = null,
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
                    const y = (i + 1) * 2 + i * H_FAC;
                    const label = item.label === "" ? "none" : item.label;
                    const isHovered =
                        hoveredKey !== null && hoveredKey === item.key;
                    const isActive =
                        activeKey !== null && activeKey === item.key;
                    const focusedKey = hoveredKey ?? activeKey;
                    const dimmed =
                        focusedKey !== null && focusedKey !== item.key;
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
                            <text
                                y={y + 9}
                                x={14}
                                aria-label={label}
                                className="text-xs fill-current"
                                style={{
                                    fontWeight: isActive ? 700 : undefined,
                                }}
                            >
                                {label}
                            </text>
                        </g>
                    );
                })}
            </g>
        </svg>
    );
}
