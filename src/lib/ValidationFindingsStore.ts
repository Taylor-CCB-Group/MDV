import { makeAutoObservable, observable, type ObservableMap } from "mobx";

export type ValidationContext = "datasource" | "chart";

export interface ValidationIssue {
    path: string;
    message: string;
    code?: string;
}

export interface ValidationFindingDetails {
    issues: ValidationIssue[];
    rawConfig?: unknown;
    details?: Record<string, unknown>;
}

export interface ValidationFinding {
    context: ValidationContext;
    /** Chart type for `chart`, datasource name for `datasource` */
    subject: string;
    signature: string;
    count: number;
    firstSeenAt: number;
    lastSeenAt: number;
    latest?: ValidationFindingDetails;
}

function issueSignature(issues: ValidationIssue[]): string {
    // Stable-ish key for grouping repeated reports of the same underlying problem.
    // Prefer code+path+message where available. Join multiple issues so zod-style multi-issue
    // errors group together.
    return issues
        .map((i) => `${i.code ?? ""}|${i.path}|${i.message}`)
        .join("\n");
}

type FindingMap = ObservableMap<string, ObservableMap<string, ValidationFinding>>;

export default class ValidationFindingsStore {
    private chartsByType: FindingMap = observable.map();
    private datasourcesByName: FindingMap = observable.map();

    constructor() {
        makeAutoObservable(this, {}, { autoBind: true });
    }

    clearAll() {
        this.chartsByType.clear();
        this.datasourcesByName.clear();
    }

    clearCharts() {
        this.chartsByType.clear();
    }

    clearDatasources() {
        this.datasourcesByName.clear();
    }

    get hasAny(): boolean {
        return this.chartsByType.size > 0 || this.datasourcesByName.size > 0;
    }

    get totalCount(): number {
        let total = 0;
        for (const [, bySig] of this.chartsByType) {
            for (const [, f] of bySig) total += f.count;
        }
        for (const [, bySig] of this.datasourcesByName) {
            for (const [, f] of bySig) total += f.count;
        }
        return total;
    }

    get chartTypeCounts(): Record<string, number> {
        const out: Record<string, number> = {};
        for (const [chartType, bySig] of this.chartsByType) {
            let count = 0;
            for (const [, f] of bySig) count += f.count;
            out[chartType] = count;
        }
        return out;
    }

    /**
     * UI/reporting-friendly structure.
     * Returned objects are derived from observable state but are plain serialisable data.
     */
    get snapshot(): {
        charts: Record<string, ValidationFinding[]>;
        datasources: Record<string, ValidationFinding[]>;
        totalCount: number;
    } {
        const charts: Record<string, ValidationFinding[]> = {};
        for (const [chartType, bySig] of this.chartsByType) {
            charts[chartType] = [...bySig.values()].sort((a, b) => b.lastSeenAt - a.lastSeenAt);
        }
        const datasources: Record<string, ValidationFinding[]> = {};
        for (const [dsName, bySig] of this.datasourcesByName) {
            datasources[dsName] = [...bySig.values()].sort((a, b) => b.lastSeenAt - a.lastSeenAt);
        }
        return { charts, datasources, totalCount: this.totalCount };
    }

    addChartFinding(chartType: string, details: ValidationFindingDetails) {
        this.addFinding("chart", chartType, details);
    }

    addDatasourceFinding(dsName: string, details: ValidationFindingDetails) {
        this.addFinding("datasource", dsName, details);
    }

    private addFinding(context: ValidationContext, subject: string, details: ValidationFindingDetails) {
        const now = Date.now();
        const signature = issueSignature(details.issues);
        const root = context === "chart" ? this.chartsByType : this.datasourcesByName;

        let bySig = root.get(subject);
        if (!bySig) {
            bySig = observable.map();
            root.set(subject, bySig);
        }

        const existing = bySig.get(signature);
        if (existing) {
            existing.count += 1;
            existing.lastSeenAt = now;
            existing.latest = details;
            return;
        }

        bySig.set(signature, {
            context,
            subject,
            signature,
            count: 1,
            firstSeenAt: now,
            lastSeenAt: now,
            latest: details,
        });
    }
}

