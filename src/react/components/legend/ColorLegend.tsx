import {
    useCallback,
    useEffect,
    useLayoutEffect,
    useRef,
    useState,
} from "react";
import type { ColorLegendSpec } from "@/react/legend/color_legend/types";
import LegendCategoricalSvg, {
    getCategoricalLabelMaxWidth,
    LEGEND_CATEGORICAL_LABEL_MAX_WIDTH,
    LEGEND_CATEGORICAL_WIDTH,
    legendCategoricalContainerHeight,
} from "@/react/components/legend/LegendCategoricalSvg";
import LegendContinuousSvg, {
    DEFAULT_CONTINUOUS_LEGEND_WIDTH,
    DEFAULT_CONTINUOUS_LEGEND_HEIGHT,
} from "@/react/components/legend/LegendContinuousSvg";

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
    const titleRef = useRef<HTMLDivElement>(null);
    const bodyRef = useRef<HTMLDivElement>(null);
    const [showTitleTooltip, setShowTitleTooltip] = useState(false);
    const [labelMaxWidth, setLabelMaxWidth] = useState(
        LEGEND_CATEGORICAL_LABEL_MAX_WIDTH,
    );
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
            const width = body.clientWidth || LEGEND_CATEGORICAL_WIDTH;
            setLabelMaxWidth(getCategoricalLabelMaxWidth(width));
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
        <div className="legend-drag-handle h-full w-full overflow-hidden">
            <div
                ref={titleRef}
                className="legend-title legend-drag-handle overflow-hidden text-ellipsis font-medium"
                title={showTitleTooltip ? label : undefined}
                style={{
                    height: "20px",
                    whiteSpace: "nowrap"
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
        const width = spec.width ?? DEFAULT_CONTINUOUS_LEGEND_WIDTH;
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
