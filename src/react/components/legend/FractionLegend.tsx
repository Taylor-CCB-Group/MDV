import { useLayoutEffect, useRef } from "react";
import type { LegendWrapperComponentProps } from "@/react/components/legend/LegendWrapper";
import type { FractionLegendSpec } from "@/react/legend/fraction_legend/types";

const FRACTION_LEGEND_WIDTH = 100;
const FRACTION_LEGEND_TITLE_HEIGHT = 20;
const FRACTION_LEGEND_ITEM_TOP = 8;
const FRACTION_LEGEND_CIRCLE_GAP = 6;
const FRACTION_LEGEND_CIRCLE_X = 40;
const FRACTION_LEGEND_TEXT_X = 70;
const FRACTION_LEGEND_VERTICAL_PADDING = 15;

function getItemPositions(items: FractionLegendSpec["items"]) {
    const positions: number[] = [];
    let contentBottom = 0;
    let currentCy = FRACTION_LEGEND_ITEM_TOP + 6;

    for (let i = 0; i < items.length; i++) {
        const radius = Number(items[i].radius) || 0;
        if (i > 0) {
            const previousRadius = Number(items[i - 1].radius) || 0;
            currentCy += previousRadius + radius + FRACTION_LEGEND_CIRCLE_GAP;
        }
        positions.push(currentCy);
        contentBottom = Math.max(contentBottom, currentCy + radius);
    }

    return {
        positions,
        bodyHeight: Math.max(
            items.length * 24 + (items.length + 2) * 2 + FRACTION_LEGEND_ITEM_TOP,
            contentBottom + 8,
        ),
    };
}

export default function FractionLegend({
    spec,
    onLayoutReady,
}: LegendWrapperComponentProps<FractionLegendSpec>) {
    const containerRef = useRef<HTMLDivElement>(null);
    const { positions, bodyHeight } = getItemPositions(spec.items);
    const autoWidth = spec.width ?? FRACTION_LEGEND_WIDTH;
    const autoHeight =
        bodyHeight +
        FRACTION_LEGEND_TITLE_HEIGHT +
        FRACTION_LEGEND_VERTICAL_PADDING;
    const width = autoWidth;
    const height = Math.min(autoHeight, spec.maxHeight ?? autoHeight);

    useLayoutEffect(() => {
        const el = containerRef.current?.parentElement;
        if (!el) {
            return;
        }
        el.style.width = `${width}px`;
        el.style.height = `${height}px`;
    }, [height, width]);

    useLayoutEffect(() => {
        onLayoutReady?.();
    }, [onLayoutReady]);

    return (
        <div ref={containerRef} className="h-full w-full overflow-hidden">
            <div
                className="legend-title legend-drag-handle overflow-hidden text-ellipsis font-medium"
                style={{
                    height: `${FRACTION_LEGEND_TITLE_HEIGHT}px`,
                    whiteSpace: "nowrap",
                }}
                title={spec.label}
            >
                {spec.label}
            </div>
            <div
                className="legend-body overflow-y-auto overflow-x-hidden w-full"
                style={{
                    height: `calc(100% - ${FRACTION_LEGEND_TITLE_HEIGHT}px)`,
                }}
            >
                <svg height={bodyHeight} width="100%" className="relative">
                    <g>
                        {spec.items.map((item, i) => {
                            const cy = positions[i];
                            return (
                                <g key={item.key}>
                                    <circle
                                        cy={cy}
                                        cx={FRACTION_LEGEND_CIRCLE_X}
                                        r={item.radius}
                                        fill="gray"
                                    />
                                    <text
                                        y={cy + 5}
                                        x={FRACTION_LEGEND_TEXT_X}
                                        alignmentBaseline="middle"
                                        className="fill-current text-xs"
                                    >
                                        {item.label}
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
