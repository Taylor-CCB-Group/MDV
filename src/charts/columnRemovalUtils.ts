import type { BaseConfig } from "./BaseChart";
import type { ChartType } from "./ChartTypes";

type ChartTypeMetadata = Pick<ChartType<BaseConfig>, "configEntriesUsingColumns"> & {
    name?: string;
};

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

function isRecord(value: unknown): value is Record<string, unknown> {
    return typeof value === "object" && value !== null;
}

export function getReferencedFields(value: unknown): string[] {
    if (value == null) {
        return [];
    }
    if (Array.isArray(value)) {
        return value.flatMap((entry) => getReferencedFields(entry));
    }
    if (typeof value === "string") {
        return [value];
    }
    if (!isRecord(value)) {
        return [];
    }
    if (typeof value.field === "string") {
        return [value.field];
    }
    if (Array.isArray(value.fields)) {
        return value.fields.filter((field): field is string => typeof field === "string");
    }
    if (typeof value.columnId === "string") {
        return [value.columnId];
    }
    if ("column" in value) {
        return getReferencedFields(value.column);
    }
    return [];
}

function getChartTitle(config: BaseConfig) {
    if (Array.isArray(config.title)) {
        return config.title.join(" ");
    }
    return config.title || config.type;
}

function getSettingUsagePaths(config: BaseConfig, chartType: ChartTypeMetadata | undefined, columnName: string) {
    const usagePaths: ColumnRemovalUsagePath[] = [];

    if (getReferencedFields((config as Record<string, unknown>).color_by).includes(columnName)) {
        usagePaths.push("color_by");
    }
    if (getReferencedFields((config as Record<string, unknown>).tooltip).includes(columnName)) {
        usagePaths.push("tooltip");
    }
    if (getReferencedFields((config as Record<string, unknown>).background_filter).includes(columnName)) {
        usagePaths.push("background_filter");
    }

    for (const entry of chartType?.configEntriesUsingColumns ?? []) {
        const value = (config as Record<string, unknown>)[entry];
        if (getReferencedFields(value).includes(columnName)) {
            usagePaths.push(entry);
        }
    }

    return usagePaths;
}

export function analyzeChartColumnImpact(
    config: BaseConfig,
    chartType: ChartTypeMetadata | undefined,
    columnName: string,
): ChartColumnImpact | null {
    const usagePaths: ColumnRemovalUsagePath[] = [];
    const paramFields = getReferencedFields(config.param);
    if (paramFields.includes(columnName)) {
        usagePaths.push("param");
    }

    const settingUsagePaths = getSettingUsagePaths(config, chartType, columnName);
    usagePaths.push(...settingUsagePaths);

    if (usagePaths.includes("param")) {
        return {
            chartId: config.id,
            chartTitle: getChartTitle(config),
            chartType: config.type,
            chartTypeLabel: chartType?.name ?? config.type,
            action: "remove_chart",
            usage: "param",
            usagePaths,
        };
    }

    if (settingUsagePaths.length === 0) {
        return null;
    }

    return {
        chartId: config.id,
        chartTitle: getChartTitle(config),
        chartType: config.type,
        chartTypeLabel: chartType?.name ?? config.type,
        action: "block_delete",
        usage: "settings",
        usagePaths: settingUsagePaths,
    };
}
