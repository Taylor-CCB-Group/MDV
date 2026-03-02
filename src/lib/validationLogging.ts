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

export function logValidationError(options: LogValidationErrorOptions): void {
    const { context, name, rawConfig, error, details } = options;
    const issues = extractIssues(error);
    const bucketName = context === "datasource" ? "datasources" : "charts";
    const label = name ?? (context === "datasource" ? "(unknown datasource)" : "(unknown chart)");

    console.error(
        `Invalid ${context} config for '${label}'. See window.mdv.validationErrors.${bucketName} for full details.`,
        issues,
    );

    if (typeof window !== "undefined" && (window as any).mdv) {
        const mdv: any = (window as any).mdv;
        const store = (mdv.validationErrors ??= { datasources: [], charts: [] });
        const bucket = store[bucketName] as any[];
        bucket.push({
            ...(name !== undefined && { name: label }),
            issues,
            rawConfig,
            ...(details && Object.keys(details).length > 0 && { details }),
        });
    }
}

