export const REGION_ID_FIELD = "region_id";
export const N_CELLS_FIELD = "n_cells";
export const N_TRANSCRIPTS_FIELD = "n_transcripts";
/** Cells-table column summed into `n_transcripts`. Omitted from the summary when the cells table has no such column. */
export const TRANSCRIPT_COUNTS_FIELD = "transcript_counts";

const RESERVED = new Set([REGION_ID_FIELD, N_CELLS_FIELD, N_TRANSCRIPTS_FIELD]);

/** `samples` unless that name is already a user table. An existing summary is detected by `n_cells`. */
export function summaryDatasourceName(existingNames: readonly string[], samplesIsSummary: boolean): string {
    if (samplesIsSummary) return "samples";
    if (existingNames.includes("sample_summary")) return "sample_summary";
    if (existingNames.includes("samples")) return "sample_summary";
    return "samples";
}

/** Text column whose category values cover the most region ids. */
export function bestIdColumn(
    regionIds: readonly string[],
    columns: readonly { field: string; values?: readonly string[] }[],
): string | null {
    const ids = new Set(regionIds);
    let best: { field: string; score: number } | null = null;
    for (const column of columns) {
        if (!column.values) continue;
        let score = 0;
        for (const value of column.values) if (ids.has(value)) score += 1;
        if (score > 0 && (!best || score > best.score)) best = { field: column.field, score };
    }
    return best?.field ?? null;
}

export function summaryFieldName(field: string): string {
    return RESERVED.has(field) ? `sample_${field}` : field;
}

export function csvEscape(value: string): string {
    if (/[",\n\r]/.test(value)) return `"${value.replace(/"/g, '""')}"`;
    return value;
}

export function buildSampleCsv(rows: readonly Record<string, string | number>[]): string {
    if (rows.length === 0) return "";
    const headers: string[] = [];
    for (const row of rows) {
        for (const key of Object.keys(row)) if (!headers.includes(key)) headers.push(key);
    }
    const lines = [headers.join(",")];
    for (const row of rows) {
        lines.push(headers.map((header) => csvEscape(String(row[header] ?? ""))).join(","));
    }
    return lines.join("\n");
}
