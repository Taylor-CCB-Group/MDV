import type { BaseConfig } from "./BaseChart";
import { chartTypes, type ChartType } from "./ChartTypes";
import type {
    AnalyzeColumnRemovalArgs,
    ChartColumnImpact,
    ColumnRemovalImpact,
    ColumnRemovalUsagePath,
    SavedViewColumnImpact,
} from "./types/columnRemoval";

type ChartTypeMetadata = Pick<ChartType<BaseConfig>, "configEntriesUsingColumns"> & {
    name?: string;
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
    const referencedFields: string[] = [];

    if (typeof value.field === "string") {
        referencedFields.push(value.field);
    } else if (isRecord(value.field)) {
        referencedFields.push(...getReferencedFields(value.field));
    }

    if (Array.isArray(value.fields)) {
        referencedFields.push(
            ...value.fields.flatMap((field) => {
                if (typeof field === "string") {
                    return [field];
                }
                if (isRecord(field)) {
                    return getReferencedFields(field);
                }
                return [];
            }),
        );
    }

    if (typeof value.columnId === "string") {
        referencedFields.push(value.columnId);
    }

    if ("column" in value) {
        referencedFields.push(...getReferencedFields(value.column));
    }

    return referencedFields;
}

function getChartTitle(config: BaseConfig) {
    if (Array.isArray(config.title)) {
        return config.title.join(" ");
    }
    return config.title || config.type;
}

function getSettingUsagePaths(config: BaseConfig, chartType: ChartTypeMetadata | undefined, columnName: string) {
    const usagePaths: ColumnRemovalUsagePath[] = [];
    const check = (path: ColumnRemovalUsagePath, value: unknown) => {
        if (getReferencedFields(value).includes(columnName)) {
            usagePaths.push(path);
        }
    };

    check("color_by", Reflect.get(config, "color_by"));
    check("tooltip", Reflect.get(config, "tooltip"));
    check("background_filter", Reflect.get(config, "background_filter"));

    for (const entry of chartType?.configEntriesUsingColumns ?? []) {
        check(entry, Reflect.get(config, entry));
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

export async function analyzeColumnRemoval({
    dataSourceName,
    columnName,
    sourceChartId = null,
    currentViewName = null,
    allViewNames = [],
    currentCharts,
    viewLoader,
}: AnalyzeColumnRemovalArgs): Promise<ColumnRemovalImpact> {
    const currentViewCharts = currentCharts
        .filter((chart) => chart.dataSourceName === dataSourceName)
        .filter((chart) => chart.config.id !== sourceChartId)
        .map((chart) =>
            analyzeChartColumnImpact(
                chart.config,
                chartTypes[chart.config.type],
                columnName,
            ),
        )
        .filter((impact) => impact !== null);

    const savedViews: SavedViewColumnImpact[] = [];
    if (viewLoader) {
        const otherViews = allViewNames.filter((viewName) => viewName !== currentViewName);
        const results = await Promise.allSettled(
            otherViews.map(async (viewName) => {
                const viewData = await viewLoader(viewName);
                const chartConfigs = viewData?.initialCharts?.[dataSourceName] ?? [];
                const charts = chartConfigs
                    .map((config) =>
                        analyzeChartColumnImpact(
                            config,
                            chartTypes[config.type],
                            columnName,
                        ),
                    )
                    .filter((impact) => impact !== null);

                if (charts.length === 0) {
                    return null;
                }

                return { viewName, charts };
            }),
        );

        const failedViews: string[] = [];
        for (const [index, result] of results.entries()) {
            if (result.status === "fulfilled" && result.value) {
                savedViews.push(result.value);
                continue;
            }
            if (result.status === "rejected") {
                failedViews.push(otherViews[index]);
            }
        }

        if (failedViews.length > 0) {
            throw new Error(
                `Failed to check column usage in saved views: ${failedViews.join(", ")}`,
            );
        }
    }

    return {
        dataSourceName,
        columnName,
        currentViewCharts,
        savedViews,
    };
}
