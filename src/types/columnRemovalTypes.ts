import type { FieldSpec, FieldSpecs } from "@/lib/columnTypeHelpers";

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
    configEntryUpdates?: Record<string, unknown>;
    tooltipUpdate?: TooltipRemovalUpdate;
    clearColorBy: boolean;
    clearBackgroundFilter: boolean;
};

export type ColumnRemovalImpact = {
    dataSource: string;
    column: string;
    charts: ChartColumnImpact[];
};

export type ParamLayoutSlot = {
    mode: "single" | "multi";
    paramIndex: number;
    entries: FieldSpecs;
};

export type TooltipConfigLike = {
    tooltip?: {
        show?: boolean;
        column?: FieldSpec | FieldSpecs;
    };
};

export type BackgroundFilterLike = {
    background_filter?: {
        column?: string;
    };
};

export type RowSummaryConfigLike = {
    type?: string;
    image?: {
        param?: number;
    };
    param?: FieldSpecs;
};

export type LegacyColorByLike = {
    color_by?: string | { column?: { field?: string } };
};
