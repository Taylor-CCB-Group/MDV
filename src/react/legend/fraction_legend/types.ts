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
};
