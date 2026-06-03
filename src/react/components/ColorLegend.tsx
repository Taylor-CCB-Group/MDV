import { useEffect, useLayoutEffect, useRef, useState } from "react";
import type { ColorLegendSpec } from "@/react/colorLegend/types";
import LegendCategoricalSvg, {
    LEGEND_CATEGORICAL_WIDTH,
    legendCategoricalContainerHeight,
} from "./LegendCategoricalSvg";
import LegendContinuousSvg, {
    DEFAULT_CONTINUOUS_LEGEND_HEIGHT,
} from "./LegendContinuousSvg";

export type ColorLegendProps = {
    spec: ColorLegendSpec;
    /** Called after DOM for this legend is committed (used to attach drag/resize on the wrapper). */
    onLayoutReady?: () => void;
    activeCategoricalValue?: string | null;
    onCategoricalItemClick?: (value: string) => void;
};

function ColorLegendCategorical({
    label,
    items,
    activeValue = null,
    onItemClick,
}: {
    label: string;
    items: { color: string; name: string; value: string }[];
    activeValue?: string | null;
    onItemClick?: (value: string) => void;
}) {
    const [hoveredValue, setHoveredValue] = useState<string | null>(null);
    const svgItems = items.map((item) => ({
        key: item.value,
        color: item.color,
        label: item.name,
    }));

    return (
        <div className="legend-drag-handle h-full w-full overflow-hidden">
            <div
                className="legend-title overflow-hidden text-ellipsis"
                title={label}
                style={{
                    height: "20px",
                    whiteSpace: "nowrap",
                }}
            >
                {label}
            </div>
            <div
                className="legend-body overflow-y-auto overflow-x-hidden w-full"
                style={{
                    height: "calc(100% - 20px)",
                }}
            >
                <LegendCategoricalSvg
                    items={svgItems}
                    hoveredKey={hoveredValue}
                    activeKey={activeValue}
                    onItemHover={setHoveredValue}
                    onItemClick={onItemClick}
                />
            </div>
        </div>
    );
}

/**
 * Inner content for the chart color legend wrapper (categorical list or continuous color bar).
 */
export default function ColorLegend({
    spec,
    onLayoutReady,
    activeCategoricalValue = null,
    onCategoricalItemClick,
}: ColorLegendProps) {
    const containerRef = useRef<HTMLDivElement>(null);

    useLayoutEffect(() => {
        spec;
        onLayoutReady?.();
    }, [spec, onLayoutReady]);

    useEffect(() => {
        const el = containerRef.current?.parentElement;
        if (!el || spec.kind !== "categorical") {
            return;
        }
        const height = legendCategoricalContainerHeight(spec.items.length);
        el.style.width = `${LEGEND_CATEGORICAL_WIDTH}px`;
        el.style.height = `${height}px`;
        el.style.border = "0.5px solid black";
    }, [spec]);

    useEffect(() => {
        const el = containerRef.current?.parentElement;
        if (!el || spec.kind !== "continuous") {
            return;
        }
        const width = spec.width ?? 120;
        const height = spec.height ?? DEFAULT_CONTINUOUS_LEGEND_HEIGHT;
        el.style.width = `${width}px`;
        el.style.height = `${height}px`;
        el.style.border = "";
    }, [spec]);

    if (spec.kind === "categorical") {
        return (
            <div ref={containerRef} className="h-full w-full">
                <ColorLegendCategorical
                    label={spec.label}
                    items={spec.items}
                    activeValue={activeCategoricalValue}
                    onItemClick={onCategoricalItemClick}
                />
            </div>
        );
    }

    return (
        <div ref={containerRef} className="h-full w-full overflow-hidden">
            <LegendContinuousSvg
                label={spec.label}
                colors={spec.colors}
                range={spec.range}
                width={spec.width}
                height={spec.height}
            />
        </div>
    );
}
