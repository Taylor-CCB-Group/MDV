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

const FRACTION_LEGEND_ITEM_GAP = 12;
const FRACTION_LEGEND_LABEL_GAP = 4;
const FRACTION_LEGEND_LEFT_PADDING = 6;
/** Approximate width of a `text-xs` character, for sizing a horizontal legend. */
const FRACTION_LEGEND_CHAR_WIDTH = 7;

type ItemLayout = { cx: number; cy: number; textX: number };

function getVerticalLayout(items: FractionLegendSpec["items"]) {
    const layout: ItemLayout[] = [];
    let contentBottom = 0;
    let currentCy = FRACTION_LEGEND_ITEM_TOP + 6;

    for (let i = 0; i < items.length; i++) {
        const radius = Number(items[i].radius) || 0;
        if (i > 0) {
            const previousRadius = Number(items[i - 1].radius) || 0;
            currentCy += previousRadius + radius + FRACTION_LEGEND_CIRCLE_GAP;
        }
        layout.push({
            cx: FRACTION_LEGEND_CIRCLE_X,
            cy: currentCy,
            textX: FRACTION_LEGEND_TEXT_X,
        });
        contentBottom = Math.max(contentBottom, currentCy + radius);
    }

    return {
        layout,
        bodyWidth: FRACTION_LEGEND_WIDTH,
        bodyHeight: Math.max(
            items.length * 24 + (items.length + 2) * 2 + FRACTION_LEGEND_ITEM_TOP,
            contentBottom + 8,
        ),
    };
}

/** Circles in one row, each followed by its label, sharing a centre line. */
function getHorizontalLayout(items: FractionLegendSpec["items"]) {
    const maxRadius = Math.max(
        0,
        ...items.map((item) => Number(item.radius) || 0),
    );
    const cy = FRACTION_LEGEND_ITEM_TOP + maxRadius;
    let x = FRACTION_LEGEND_LEFT_PADDING;
    const layout = items.map((item) => {
        const radius = Number(item.radius) || 0;
        const cx = x + radius;
        const textX = cx + radius + FRACTION_LEGEND_LABEL_GAP;
        x =
            textX +
            item.label.length * FRACTION_LEGEND_CHAR_WIDTH +
            FRACTION_LEGEND_ITEM_GAP;
        return { cx, cy, textX };
    });
    return { layout, bodyWidth: x, bodyHeight: cy + maxRadius + 6 };
}

function getItemLayout(spec: FractionLegendSpec) {
    return spec.orientation === "horizontal"
        ? getHorizontalLayout(spec.items)
        : getVerticalLayout(spec.items);
}

/** Size the legend sets on its wrapper, so charts can lay out around it before it renders. */
export function getFractionLegendSize(spec: FractionLegendSpec) {
    const { bodyWidth, bodyHeight } = getItemLayout(spec);
    const titleWidth =
        spec.label.length * FRACTION_LEGEND_CHAR_WIDTH +
        FRACTION_LEGEND_LEFT_PADDING * 2;
    const autoWidth =
        spec.orientation === "horizontal"
            ? Math.max(bodyWidth, titleWidth)
            : FRACTION_LEGEND_WIDTH;
    const autoHeight =
        bodyHeight +
        FRACTION_LEGEND_TITLE_HEIGHT +
        FRACTION_LEGEND_VERTICAL_PADDING;
    return {
        width: spec.width ?? autoWidth,
        height: Math.min(autoHeight, spec.maxHeight ?? autoHeight),
    };
}

export default function FractionLegend({
    spec,
    onLayoutReady,
}: LegendWrapperComponentProps<FractionLegendSpec>) {
    const containerRef = useRef<HTMLDivElement>(null);
    const { layout, bodyHeight } = getItemLayout(spec);
    const { width, height } = getFractionLegendSize(spec);

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
                            const { cx, cy, textX } = layout[i];
                            return (
                                <g key={item.key}>
                                    <circle
                                        cy={cy}
                                        cx={cx}
                                        r={item.radius}
                                        fill="gray"
                                    />
                                    <text
                                        y={cy + 5}
                                        x={textX}
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
