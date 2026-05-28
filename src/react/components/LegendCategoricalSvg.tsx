export type LegendCategoricalSvgItem = {
    key: string;
    color: string;
    label: string;
};

export type LegendCategoricalSvgProps = {
    items: LegendCategoricalSvgItem[];
    /** Called when pointer enters/leaves a row (optional, used by field contour legend). */
    onItemHover?: (key: string | null) => void;
    hoveredKey?: string | null;
};

const H_FAC = 12;
const MAX_HEIGHT = 215;

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
    hoveredKey = null,
}: LegendCategoricalSvgProps) {
    const height = legendCategoricalBodyHeight(items.length);

    return (
        <svg height={height} width={180} className="relative">
            <g>
                {items.map((item, i) => {
                    const y = (i + 1) * 2 + i * H_FAC;
                    const isHovered =
                        hoveredKey !== null && hoveredKey === item.key;
                    const dimmed =
                        hoveredKey !== null && hoveredKey !== item.key;
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
                            className={onItemHover ? "cursor-pointer" : undefined}
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
                                stroke={isHovered ? "currentColor" : "none"}
                                strokeWidth={isHovered ? 2 : 0}
                            />
                            <text
                                y={y + 9}
                                x={14}
                                className="text-xs fill-current"
                            >
                                {item.label === "" ? "none" : item.label}
                            </text>
                        </g>
                    );
                })}
            </g>
        </svg>
    );
}
