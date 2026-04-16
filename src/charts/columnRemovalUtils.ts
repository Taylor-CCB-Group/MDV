import type { BaseConfig } from "./BaseChart";
import type { ChartType } from "./ChartTypes";
import { flattenFields, isMultiColumn, type FieldSpec, type FieldSpecs } from "@/lib/columnTypeHelpers";

export type ChartColumnImpactAction = "update" | "delete";

export type ChartColumnImpactReasonKind =
    | "param.single"
    | "param.multi"
    | "color_by"
    | "tooltip.single"
    | "tooltip.multi"
    | "background_filter"
    | "config_entry.single"
    | "config_entry.multi";

export type ChartColumnImpactReason = {
    kind: ChartColumnImpactReasonKind;
    entry?: string;
};

export type ParamSlotImpact = {
    paramIndex: number;
    mode: "single" | "multi";
    action: "delete" | "prune";
    remainingFields: number;
};

export type TooltipRemovalUpdate = {
    nextColumn?: FieldSpec | FieldSpecs;
    disableTooltip: boolean;
};

export type ChartColumnImpact = {
    chartId: string;
    chartTitle: string;
    chartType: string;
    chartTypeLabel: string;
    isSourceChart: boolean;
    action: ChartColumnImpactAction;
    reasons: ChartColumnImpactReason[];
    paramSlotImpacts: ParamSlotImpact[];
    nextParam?: FieldSpecs;
    configEntryUpdates?: Record<string, FieldSpec | FieldSpecs | undefined>;
    tooltipUpdate?: TooltipRemovalUpdate;
    clearColorBy: boolean;
    clearBackgroundFilter: boolean;
};

export type ColumnRemovalImpact = {
    dataSource: string;
    column: string;
    charts: ChartColumnImpact[];
};

type ParamLayoutSlot = {
    mode: "single" | "multi";
    paramIndex: number;
    entries: FieldSpecs;
};

type TooltipConfigLike = {
    tooltip?: {
        show?: boolean;
        column?: FieldSpec | FieldSpecs;
    };
};

type BackgroundFilterLike = {
    background_filter?: {
        column?: string;
    };
};

type RowSummaryConfigLike = {
    type?: string;
    image?: {
        param?: number;
    };
    param?: FieldSpecs;
};

type LegacyColorByLike = {
    color_by?: string | { column?: { field?: string } };
};

function getFieldSpecs(value: unknown): FieldSpecs | null {
    if (value == null) {
        return null;
    }
    return Array.isArray(value) ? value : [value as FieldSpec];
}

function containsField(value: unknown, column: string): boolean {
    const specs = getFieldSpecs(value);
    return specs ? specs.some((spec) => flattenFields(spec).includes(column)) : false;
}

function getColorField(config: BaseConfig): string | undefined {
    const colorBy = (config as LegacyColorByLike).color_by;
    if (!colorBy) {
        return undefined;
    }
    if (typeof colorBy === "string") {
        return colorBy;
    }
    return colorBy.column?.field;
}

function getParamLayouts(
    param: FieldSpecs | undefined,
    chartType: Pick<ChartType<BaseConfig>, "params"> | undefined,
): ParamLayoutSlot[] {
    const currentParam =
        param == null ? [] : Array.isArray(param) ? param : [param];
    const chartParams = chartType?.params ?? [];

    if (currentParam.length === 0) {
        return [];
    }
    if (chartParams.length === 0) {
        return currentParam.map((entry, index) => ({
            mode: "single",
            paramIndex: index,
            entries: [entry],
        }));
    }

    const multiIndex = chartParams.findIndex((slot) => isMultiColumn(slot.type));
    if (multiIndex === -1) {
        return currentParam.map((entry, index) => ({
            mode: "single",
            paramIndex: index,
            entries: [entry],
        }));
    }

    if (multiIndex === 0) {
        const singleSuffixCount = chartParams.length - 1;
        const multiEnd = Math.max(0, currentParam.length - singleSuffixCount);
        return [
            {
                mode: "multi",
                paramIndex: 0,
                entries: currentParam.slice(0, multiEnd),
            },
            ...chartParams.slice(1).map((_, index) => ({
                mode: "single" as const,
                paramIndex: index + 1,
                entries: currentParam[multiEnd + index] ? [currentParam[multiEnd + index]] : [],
            })),
        ];
    }

    if (multiIndex === chartParams.length - 1) {
        const singlePrefixCount = chartParams.length - 1;
        return [
            ...chartParams.slice(0, singlePrefixCount).map((_, index) => ({
                mode: "single" as const,
                paramIndex: index,
                entries: currentParam[index] ? [currentParam[index]] : [],
            })),
            {
                mode: "multi",
                paramIndex: multiIndex,
                entries: currentParam.slice(singlePrefixCount),
            },
        ];
    }

    return currentParam.map((entry, index) => ({
        mode: "single",
        paramIndex: index,
        entries: [entry],
    }));
}

function getTooltipUpdate(config: BaseConfig, column: string): TooltipRemovalUpdate | undefined {
    const tooltipColumn = (config as TooltipConfigLike).tooltip?.column;
    if (!tooltipColumn || !containsField(tooltipColumn, column)) {
        return undefined;
    }

    if (Array.isArray(tooltipColumn)) {
        const nextColumn = tooltipColumn.filter((entry) => !flattenFields(entry).includes(column));
        return {
            nextColumn: nextColumn.length > 0 ? nextColumn : undefined,
            disableTooltip: nextColumn.length === 0,
        };
    }

    return {
        nextColumn: undefined,
        disableTooltip: true,
    };
}

export function analyzeChartColumnImpact(
    config: BaseConfig,
    chartType: Pick<ChartType<BaseConfig>, "params" | "configEntriesUsingColumns"> & { name?: string } | undefined,
    column: string,
    opts: { sourceChartId?: string } = {},
): ChartColumnImpact | null {
    const reasons: ChartColumnImpactReason[] = [];
    const paramSlotImpacts: ParamSlotImpact[] = [];
    const nextParamSlots: ParamLayoutSlot[] = [];
    let action: ChartColumnImpactAction = "update";
    let touched = false;

    // This is the shared prune-vs-delete decision point. Each chart's `params`
    // metadata tells us whether a removed field hits a required single slot or
    // a prunable multi-column slot.
    for (const slot of getParamLayouts(config.param, chartType)) {
        if (slot.entries.length === 0) {
            nextParamSlots.push(slot);
            continue;
        }

        const matchingEntries = slot.entries.filter((entry) => flattenFields(entry).includes(column));
        if (matchingEntries.length === 0) {
            nextParamSlots.push(slot);
            continue;
        }

        touched = true;
        if (slot.mode === "single") {
            reasons.push({ kind: "param.single" });
            paramSlotImpacts.push({
                paramIndex: slot.paramIndex,
                mode: "single",
                action: "delete",
                remainingFields: 0,
            });
            action = "delete";
            nextParamSlots.push(slot);
            continue;
        }

        const remainingEntries = slot.entries.filter((entry) => !flattenFields(entry).includes(column));
        reasons.push({ kind: "param.multi" });
        paramSlotImpacts.push({
            paramIndex: slot.paramIndex,
            mode: "multi",
            action: remainingEntries.length === 0 ? "delete" : "prune",
            remainingFields: remainingEntries.length,
        });
        if (remainingEntries.length === 0) {
            action = "delete";
        }
        nextParamSlots.push({
            ...slot,
            entries: remainingEntries,
        });
    }

    const clearColorBy = getColorField(config) === column;
    if (clearColorBy) {
        touched = true;
        reasons.push({ kind: "color_by" });
    }

    const tooltipUpdate = getTooltipUpdate(config, column);
    if (tooltipUpdate) {
        touched = true;
        reasons.push({
            kind: Array.isArray((config as TooltipConfigLike).tooltip?.column)
                ? "tooltip.multi"
                : "tooltip.single",
        });
    }

    const clearBackgroundFilter =
        (config as BackgroundFilterLike).background_filter?.column === column;
    if (clearBackgroundFilter) {
        touched = true;
        reasons.push({ kind: "background_filter" });
    }

    // `configEntriesUsingColumns` covers column references that are not part of
    // `config.param`, for example image labels or sort-by fields.
    const configEntryUpdates: Record<string, FieldSpec | FieldSpecs | undefined> = {};
    for (const entry of chartType?.configEntriesUsingColumns ?? []) {
        const value = (config as Record<string, unknown>)[entry];
        if (value == null || !containsField(value, column)) {
            continue;
        }

        touched = true;
        if (Array.isArray(value)) {
            const nextValue = value.filter((item) => !flattenFields(item as FieldSpec).includes(column));
            reasons.push({ kind: "config_entry.multi", entry });
            if (nextValue.length === 0) {
                action = "delete";
            } else {
                configEntryUpdates[entry] = nextValue as FieldSpecs;
            }
            continue;
        }

        reasons.push({ kind: "config_entry.single", entry });
        action = "delete";
    }

    if (!touched) {
        return null;
    }

    // `row_summary_box` stores an image key as an index into `param`, so
    // deleting that specific field is destructive even though the chart also
    // has a multi-column param.
    const rowSummaryConfig = config as RowSummaryConfigLike;
    if (rowSummaryConfig.type === "row_summary_box" && rowSummaryConfig.image) {
        const imageParamIndex = rowSummaryConfig.image.param;
        const currentParam = Array.isArray(rowSummaryConfig.param) ? rowSummaryConfig.param : [];
        if (
            imageParamIndex !== undefined &&
            currentParam[imageParamIndex] &&
            flattenFields(currentParam[imageParamIndex]).includes(column)
        ) {
            action = "delete";
        }
    }

    return {
        chartId: config.id,
        chartTitle: config.title,
        chartType: config.type,
        chartTypeLabel: chartType?.name ?? config.type,
        isSourceChart: config.id === opts.sourceChartId,
        action,
        reasons,
        paramSlotImpacts,
        nextParam: nextParamSlots.flatMap((slot) => slot.entries),
        configEntryUpdates: Object.keys(configEntryUpdates).length > 0 ? configEntryUpdates : undefined,
        tooltipUpdate,
        clearColorBy,
        clearBackgroundFilter,
    };
}

export function describeColumnImpactReason(reason: ChartColumnImpactReason): string {
    switch (reason.kind) {
        case "param.single":
            return "parameter";
        case "param.multi":
            return "multiple parameter";
        case "color_by":
            return "color by";
        case "tooltip.single":
        case "tooltip.multi":
            return "tooltip";
        case "background_filter":
            return "background filter";
        case "config_entry.single":
        case "config_entry.multi":
            return reason.entry ? `config: ${reason.entry}` : "config entry";
    }
}
