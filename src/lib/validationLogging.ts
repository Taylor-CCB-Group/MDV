import type ValidationFindingsStore from "./ValidationFindingsStore";
import type { ChartConfig } from "@/charts/schemas/ChartConfigSchema";

type ValidationContext = "datasource" | "chart";

export interface ValidationIssue {
    path: string;
    message: string;
    code?: string;
}

export interface LogValidationErrorOptions {
    context: ValidationContext;
    name?: string;
    rawConfig?: unknown;
    error: unknown;
    details?: Record<string, unknown>;
}

function extractIssues(error: unknown): ValidationIssue[] {
    const err: any = error;
    const entries: any[] | null =
        Array.isArray(err?.issues) ? err.issues :
        Array.isArray(err?.errors) ? err.errors :
        null;

    if (entries) {
        return entries.map((issue) => ({
            path: Array.isArray(issue.path) ? issue.path.join(".") : "",
            message: issue.message ?? String(issue),
            code: issue.code,
        }));
    }

    const message =
        error instanceof Error ? error.message :
        typeof error === "string" ? error :
        JSON.stringify(error);

    return [{
        path: "",
        message,
    }];
}

function getFindingsStore(): ValidationFindingsStore | undefined {
    if (typeof window === "undefined") return undefined;
    const mdv = window.mdv;
    if (!mdv) return undefined;
    const cm = mdv.chartManager;
    if (!cm) return undefined;
    return cm.validationFindings;
}

export function logValidationError(options: LogValidationErrorOptions) {
    const { context, name, rawConfig, error, details } = options;
    const issues = extractIssues(error);
    const bucketName = context === "datasource" ? "datasources" : "charts";
    const label = name ?? (context === "datasource" ? "(unknown datasource)" : "(unknown chart)");

    console.error(
        `Invalid ${context} config for '${label}'. See Debug / Report for aggregated details.`,
        issues,
    );

    const store = getFindingsStore();
    if (!store) return;

    if (context === "chart") {
        store.addChartFinding(label, {
            issues,
            rawConfig,
            ...(details && Object.keys(details).length > 0 ? { details } : {}),
        });
        return;
    }

    store.addDatasourceFinding(label, {
        issues,
        rawConfig,
        ...(details && Object.keys(details).length > 0 ? { details } : {}),
    });
}

export function logChartValidationError(
    config: ChartConfig,
    error: unknown,
    options?: { details?: Record<string, unknown> },
) {
    const cfg = config as { type?: string; version?: string };
    const label = cfg.version ? `${cfg.type ?? "unknown"} v${cfg.version}` : (cfg.type ?? "unknown");
    logValidationError({
        context: "chart",
        name: label,
        rawConfig: config,
        error,
        details: options?.details,
    });
}

export function useValidationErrors() {

}