import type { BaseConfig } from "../BaseChart";

export type ColumnRemovalAction = "remove_chart" | "block_delete";
export type ColumnRemovalUsage = "param" | "settings";
export type ColumnRemovalUsagePath =
    | "param"
    | "color_by"
    | "tooltip"
    | "background_filter"
    | string;

export type ChartColumnImpact = {
    chartId?: string;
    chartTitle: string;
    chartType: string;
    chartTypeLabel: string;
    action: ColumnRemovalAction;
    usage: ColumnRemovalUsage;
    usagePaths: ColumnRemovalUsagePath[];
};

export type SavedViewColumnImpact = {
    viewName: string;
    charts: ChartColumnImpact[];
};

export type ColumnRemovalImpact = {
    dataSourceName: string;
    columnName: string;
    currentViewCharts: ChartColumnImpact[];
    savedViews: SavedViewColumnImpact[];
};

export type CurrentChartInfo = {
    dataSourceName: string;
    config: BaseConfig;
};

type LoadedViewData = {
    initialCharts?: Record<string, BaseConfig[]>;
} | null | undefined;

export type ViewLoader = (viewName: string) => Promise<LoadedViewData>;

export type AnalyzeColumnRemovalArgs = {
    dataSourceName: string;
    columnName: string;
    sourceChartId?: string | null;
    currentViewName?: string | null;
    allViewNames?: string[];
    currentCharts: CurrentChartInfo[];
    viewLoader?: ViewLoader;
};
