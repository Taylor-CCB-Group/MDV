import type {
    BackgroundFilterLike,
    ChartColumnImpact,
    ChartColumnImpactAction,
    ChartColumnImpactReason,
    LegacyColorByLike,
    ParamLayoutSlot,
    ParamSlotImpact,
    RowSummaryConfigLike,
    TooltipConfigLike,
    TooltipRemovalUpdate,
} from "@/types/columnRemovalTypes";
import type { BaseConfig } from "./BaseChart";
import type { ChartType } from "./ChartTypes";
import { flattenFields, isMultiColumn, type FieldSpec, type FieldSpecs } from "@/lib/columnTypeHelpers";

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

function getColorFields(config: BaseConfig): string[] {
    const colorBy = (config as LegacyColorByLike).color_by;
    if (!colorBy) {
        return [];
    }
    if (typeof colorBy === "string") {
        return [colorBy];
    }
    if ("column" in colorBy) {
        return colorBy.column?.field ? [colorBy.column.field] : [];
    }
    return flattenFields(colorBy as FieldSpec);
}

function getParamLayouts(
    param: FieldSpecs | undefined,
    chartType: (Pick<ChartType<BaseConfig>, "params"> & { name?: string }) | undefined,
): ParamLayoutSlot[] {
    const currentParam = param == null ? [] : Array.isArray(param) ? param : [param];
    const chartParams = chartType?.params ?? [];

    if (currentParam.length === 0) {
        return [];
    }
    if (chartParams.length === 0) {
        // No metadata means we cannot infer any optional multi-column slot, so
        // the safest assumption is that every param entry is required.
        return currentParam.map((entry, index) => ({
            mode: "single",
            paramIndex: index,
            entries: [entry],
        }));
    }

    const multiIndex = chartParams.findIndex((slot) => isMultiColumn(slot.type));
    if (multiIndex === -1) {
        // Simple case: config.param is just a positional list of required
        // single-column params.
        return currentParam.map((entry, index) => ({
            mode: "single",
            paramIndex: index,
            entries: [entry],
        }));
    }

    if (multiIndex === 0) {
        const singleSuffixCount = chartParams.length - 1;
        const multiEnd = Math.max(0, currentParam.length - singleSuffixCount);
        // Example layout:
        //   params metadata: [multi, single, single]
        //   config.param:     [a, b, c, d]
        // becomes:
        //   slot0(multi): [a, b]
        //   slot1(single): [c]
        //   slot2(single): [d]
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
        // Example layout:
        //   params metadata: [single, multi]
        //   config.param:     [a, b, c]
        // becomes:
        //   slot0(single): [a]
        //   slot1(multi):  [b, c]
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

    // The UI and shared cleanup code only understand one multi-column slot,
    // and only when it is first or last. Throw here rather than silently
    // guessing a different layout and pruning the wrong field.
    throw new Error(
        `Unsupported column-removal param layout for '${chartType?.name ?? "unknown chart"}': multi-column params must be first or last.`,
    );
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
    chartType: (Pick<ChartType<BaseConfig>, "params" | "configEntriesUsingColumns"> & { name?: string }) | undefined,
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
    //
    // Once `action` flips to `"delete"` we keep it that way. Later config
    // references may still add useful reason metadata, but they do not make the
    // chart valid again.
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

    const clearColorBy = getColorFields(config).includes(column);
    if (clearColorBy) {
        // Secondary color binding is survivable: clear it and let the chart
        // fall back to its default coloring behavior.
        touched = true;
        reasons.push({ kind: "color_by" });
    }

    const tooltipUpdate = getTooltipUpdate(config, column);
    if (tooltipUpdate) {
        touched = true;
        reasons.push({
            kind: Array.isArray((config as TooltipConfigLike).tooltip?.column) ? "tooltip.multi" : "tooltip.single",
        });
    }

    const clearBackgroundFilter = (config as BackgroundFilterLike).background_filter?.column === column;
    if (clearBackgroundFilter) {
        // Background filter is also a survivable secondary reference.
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
        // Exclude unaffected charts from the impact list entirely so both the
        // dialog and execution path stay focused only on real dependents.
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
        // Convert the surviving logical slots back into the flat param array
        // shape expected by `config.param` / `setParams()`.
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
