import {
    CONTINUOUS_HORIZONTAL_PADDING,
    CONTINUOUS_MIN_LAYOUT_WIDTH,
    CONTINUOUS_TICK_TARGET_SPACING,
    LEGEND_CATEGORICAL_LABEL_PADDING_X,
    LEGEND_CATEGORICAL_LABEL_START_X,
    LEGEND_CATEGORICAL_MAX_BODY_HEIGHT,
    LEGEND_CATEGORICAL_ROW_HEIGHT,
    LEGEND_CATEGORICAL_WIDTH,
    LEGEND_LABEL_FONT,
} from "./legendConstants";
import type { ContinuousLegendLayout, GradientStop } from "./legendTypes";

const ellipsis = "…";
let measureCanvas: HTMLCanvasElement | null = null;

function getMeasureContext(): CanvasRenderingContext2D | null {
    if (typeof document === "undefined") {
        return null;
    }
    if (
        typeof navigator !== "undefined" &&
        navigator.userAgent.includes("jsdom")
    ) {
        return null;
    }
    measureCanvas ??= document.createElement("canvas");
    try {
        return measureCanvas.getContext("2d");
    } catch {
        return null;
    }
}

export function measureLegendLabelWidth(
    text: string,
    font = LEGEND_LABEL_FONT,
): number {
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
    while (
        display.length > 0 &&
        measureLegendLabelWidth(`${display}${ellipsis}`) > maxWidth
    ) {
        display = display.slice(0, -1);
    }
    return {
        display: `${display}${ellipsis}`,
        truncated: true,
        full,
    };
}

export function formatContinuousTick(value: unknown): string {
    const num = Number(value);
    if (!Number.isFinite(num)) {
        return String(value);
    }
    if (Math.abs(num) >= 10000 || (Math.abs(num) > 0 && Math.abs(num) < 0.01)) {
        return num.toExponential(1).replace("+", "");
    }
    return String(num);
}

export const LEGEND_CATEGORICAL_LABEL_MAX_WIDTH =
    getCategoricalLabelMaxWidth(LEGEND_CATEGORICAL_WIDTH);

export function getCategoricalLabelMaxWidth(width: number): number {
    return Math.max(
        0,
        width -
            LEGEND_CATEGORICAL_LABEL_START_X -
            LEGEND_CATEGORICAL_LABEL_PADDING_X,
    );
}

export function legendCategoricalRowY(index: number): number {
    return (index + 1) * 2 + index * LEGEND_CATEGORICAL_ROW_HEIGHT;
}

export function legendCategoricalBodyHeight(itemCount: number): number {
    return (
        itemCount * LEGEND_CATEGORICAL_ROW_HEIGHT + (itemCount + 2) * 2
    );
}

export function legendCategoricalContainerHeight(itemCount: number): number {
    const bodyHeight = legendCategoricalBodyHeight(itemCount);
    return bodyHeight > LEGEND_CATEGORICAL_MAX_BODY_HEIGHT
        ? LEGEND_CATEGORICAL_MAX_BODY_HEIGHT + 35
        : bodyHeight + 35;
}

export function getContinuousLegendLayout(
    width: number,
    hasLabel: boolean,
): ContinuousLegendLayout {
    const layoutWidth = Math.max(width, CONTINUOUS_MIN_LAYOUT_WIDTH);
    const axisWidth = layoutWidth - CONTINUOUS_HORIZONTAL_PADDING * 2;

    return {
        layoutWidth,
        barX: CONTINUOUS_HORIZONTAL_PADDING,
        barY: hasLabel ? 15 : 0,
        barHeight: 10,
        axisWidth,
        labelMaxWidth: Math.max(0, axisWidth),
        tickCount: Math.max(
            2,
            Math.ceil(layoutWidth / CONTINUOUS_TICK_TARGET_SPACING),
        ),
    };
}

export function getGradientStops(colors: string[]): GradientStop[] {
    return colors.map((color, index) => {
        const pct =
            colors.length <= 1 ? 0 : (index / (colors.length - 1)) * 100;
        return {
            offset: `${pct}%`,
            color,
            index,
        };
    });
}
