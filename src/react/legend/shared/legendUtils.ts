import { scaleLinear } from "d3-scale";
import { dateTickFormat } from "@/lib/dateFormat";
import {
    CONTINUOUS_COMPACT_BAR_HEIGHT,
    CONTINUOUS_COMPACT_LABEL_HEIGHT,
    CONTINUOUS_COMPACT_MAX_BAR_WIDTH,
    CONTINUOUS_COMPACT_TICK_LENGTH,
    CONTINUOUS_COMPACT_TICK_SPACING,
    CONTINUOUS_COMPACT_TITLE_HEIGHT,
    CONTINUOUS_HORIZONTAL_PADDING,
    CONTINUOUS_MIN_LAYOUT_WIDTH,
    CONTINUOUS_TICK_TARGET_SPACING,
    CONTINUOUS_TITLE_GRADIENT_GAP,
    CONTINUOUS_TITLE_HEIGHT,
    DEFAULT_CONTINUOUS_LEGEND_WIDTH,
    LEGEND_CATEGORICAL_LABEL_PADDING_X,
    LEGEND_CATEGORICAL_LABEL_START_X,
    LEGEND_CATEGORICAL_MAX_BODY_HEIGHT,
    LEGEND_CATEGORICAL_ROW_HEIGHT,
    LEGEND_CATEGORICAL_WIDTH,
    LEGEND_LABEL_FONT,
} from "@/react/legend/shared/legendConstants";
import type { ContinuousLegendLayout, GradientStop } from "@/react/legend/shared/legendTypes";

const ellipsis = "…";
let measureCanvas: HTMLCanvasElement | null = null;

function getMeasureContext(): CanvasRenderingContext2D | null {
    if (typeof document === "undefined") {
        return null;
    }
    if (typeof navigator !== "undefined" && navigator.userAgent.includes("jsdom")) {
        return null;
    }
    measureCanvas ??= document.createElement("canvas");
    try {
        return measureCanvas.getContext("2d");
    } catch {
        return null;
    }
}

export function measureLegendLabelWidth(text: string, font = LEGEND_LABEL_FONT): number {
    const ctx = getMeasureContext();
    if (!ctx) {
        return text.length * 7;
    }
    ctx.font = font;
    return ctx.measureText(text).width;
}

export function formatLegendLabel(
    label: string,
    maxWidth: number,
): { display: string; truncated: boolean; full: string } {
    const full = label === "" ? "none" : label;
    if (maxWidth <= 0 || measureLegendLabelWidth(full) <= maxWidth) {
        return { display: full, truncated: false, full };
    }

    let display = full;
    while (display.length > 0 && measureLegendLabelWidth(`${display}${ellipsis}`) > maxWidth) {
        display = display.slice(0, -1);
    }
    return {
        display: `${display}${ellipsis}`,
        truncated: true,
        full,
    };
}

export function formatContinuousTick(value: unknown, isDate = false): string {
    const num = Number(value);
    if (!Number.isFinite(num)) {
        return String(value);
    }
    if (isDate) {
        return dateTickFormat(num);
    }
    if (Math.abs(num) >= 10000 || (Math.abs(num) > 0 && Math.abs(num) < 0.01)) {
        return num.toExponential(1).replace("+", "");
    }
    return String(num);
}

export const LEGEND_CATEGORICAL_LABEL_MAX_WIDTH = getCategoricalLabelMaxWidth(LEGEND_CATEGORICAL_WIDTH);

export function getCategoricalLabelMaxWidth(width: number): number {
    return Math.max(0, width - LEGEND_CATEGORICAL_LABEL_START_X - LEGEND_CATEGORICAL_LABEL_PADDING_X);
}

export function legendCategoricalRowY(index: number): number {
    return (index + 1) * 2 + index * LEGEND_CATEGORICAL_ROW_HEIGHT;
}

export function legendCategoricalBodyHeight(itemCount: number): number {
    return itemCount * LEGEND_CATEGORICAL_ROW_HEIGHT + (itemCount + 2) * 2;
}

export function legendCategoricalContainerHeight(itemCount: number): number {
    const bodyHeight = legendCategoricalBodyHeight(itemCount);
    return bodyHeight > LEGEND_CATEGORICAL_MAX_BODY_HEIGHT ? LEGEND_CATEGORICAL_MAX_BODY_HEIGHT + 35 : bodyHeight + 35;
}

export function getContinuousLegendLayout(width: number, hasLabel: boolean, compact = false): ContinuousLegendLayout {
    const layoutWidth = Math.max(width, CONTINUOUS_MIN_LAYOUT_WIDTH);
    if (compact) {
        const axisWidth = Math.min(layoutWidth - CONTINUOUS_HORIZONTAL_PADDING * 2, CONTINUOUS_COMPACT_MAX_BAR_WIDTH);
        return {
            layoutWidth,
            barX: CONTINUOUS_HORIZONTAL_PADDING,
            barY: hasLabel ? CONTINUOUS_COMPACT_TITLE_HEIGHT : 0,
            barHeight: CONTINUOUS_COMPACT_BAR_HEIGHT,
            axisWidth,
            labelMaxWidth: Math.max(0, layoutWidth - 10),
            tickCount: Math.max(2, Math.floor(axisWidth / CONTINUOUS_COMPACT_TICK_SPACING)),
            tickLength: CONTINUOUS_COMPACT_TICK_LENGTH,
        };
    }
    const axisWidth = layoutWidth - CONTINUOUS_HORIZONTAL_PADDING * 2;

    return {
        layoutWidth,
        barX: CONTINUOUS_HORIZONTAL_PADDING,
        barY: hasLabel ? CONTINUOUS_TITLE_HEIGHT + CONTINUOUS_TITLE_GRADIENT_GAP : 0,
        barHeight: 10,
        axisWidth,
        labelMaxWidth: Math.max(0, axisWidth),
        tickCount: Math.max(2, Math.ceil(layoutWidth / CONTINUOUS_TICK_TARGET_SPACING)),
        tickLength: CONTINUOUS_AXIS_TICK_LENGTH,
    };
}

const CONTINUOUS_AXIS_TICK_LENGTH = 6;
const CONTINUOUS_AXIS_LABEL_PADDING = 8;
const CONTINUOUS_LEGEND_MIN_HEIGHT = 64;

export type ContinuousLegendHeightOptions = {
    width?: number;
    hasLabel?: boolean;
    isDate?: boolean;
    /** The strip's slimmer form: tick labels are level, so the height doesn't depend on them. */
    compact?: boolean;
};

export function getContinuousLegendContainerHeight(
    range: [number, number],
    options: ContinuousLegendHeightOptions = {},
): number {
    const width = options.width ?? DEFAULT_CONTINUOUS_LEGEND_WIDTH;
    const hasLabel = options.hasLabel ?? true;
    const isDate = options.isDate ?? false;
    if (options.compact) {
        const layout = getContinuousLegendLayout(width, hasLabel, true);
        return layout.barY + layout.barHeight + layout.tickLength + CONTINUOUS_COMPACT_LABEL_HEIGHT;
    }
    const layout = getContinuousLegendLayout(width, hasLabel);

    const scale = scaleLinear().domain([range[0], range[1]]).range([0, layout.axisWidth]);
    const ticks = scale.ticks(layout.tickCount);

    let longestTickWidth = 0;
    for (const tick of ticks) {
        longestTickWidth = Math.max(longestTickWidth, measureLegendLabelWidth(formatContinuousTick(tick, isDate)));
    }

    // Matches axisBottom labels: text-anchor end, dx -.4em, dy .4em, rotate(-45).
    const fontSize = 12;
    const labelDrop = Math.ceil((longestTickWidth + fontSize) * Math.SQRT1_2 + CONTINUOUS_AXIS_LABEL_PADDING);

    const contentHeight = layout.barY + layout.barHeight + CONTINUOUS_AXIS_TICK_LENGTH + labelDrop;

    return Math.max(CONTINUOUS_LEGEND_MIN_HEIGHT, contentHeight);
}

export function getGradientStops(colors: string[]): GradientStop[] {
    return colors.map((color, index) => {
        const pct = colors.length <= 1 ? 0 : (index / (colors.length - 1)) * 100;
        return {
            offset: `${pct}%`,
            color,
            index,
        };
    });
}
