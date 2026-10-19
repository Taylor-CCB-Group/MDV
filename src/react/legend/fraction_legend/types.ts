export type FractionLegendItem = {
    key: string;
    label: string;
    radius: number;
};

export type FractionLegendSpec = {
    label: string;
    items: FractionLegendItem[];
    width?: number;
    maxHeight?: number;
    /** Circles stacked in a column (default) or side by side in a row. */
    orientation?: "vertical" | "horizontal";
};
