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
    labelMaxWidth?: number;
};

export type LegendContinuousSvgProps = {
    label: string;
    colors: string[];
    range: [number, number];
    width?: number;
    height?: number;
    activeRange?: [number, number] | null;
    onRangeChange?: (range: [number, number] | null) => void;
};

export type ContinuousLegendLayout = {
    layoutWidth: number;
    barX: number;
    barY: number;
    barHeight: number;
    axisWidth: number;
    labelMaxWidth: number;
    tickCount: number;
};

export type GradientStop = {
    offset: string;
    color: string;
    index: number;
};
