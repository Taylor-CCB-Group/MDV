/**
 * Brings the open page's datastores up to date with the columns jobs wrote on the server, so new
 * columns can be used without a reload.
 *
 * Ingest stamps each output column's metadata with a `provenance` pointer naming the job that wrote
 * it. A column the datastore lacks is added; a column whose pointer names a different job was
 * overwritten by a re-run, so its metadata is replaced and any loaded data reloaded. Running the
 * sync again changes nothing.
 */
import type DataStore from "@/datastore/DataStore";

/** A column's metadata as `/get_configs` returns it. */
export type ServerColumn = { field: string; deleted?: boolean; provenance?: { kind?: string; job_id?: string } };
export type ServerDataSource = { name: string; columns?: ServerColumn[] };

/** The parts of a ChartManager that `syncJobColumns` uses. */
export type ColumnSyncTarget = {
    dataSources: { name: string; dataStore: DataStore }[];
    loadColumnSet(columns: string[], dataSource: string, callback: () => void): void;
};

/** The job that wrote the column, from its provenance pointer; undefined for any other column. */
function jobId(column: object): string | undefined {
    const p = "provenance" in column ? column.provenance : undefined;
    if (typeof p !== "object" || p === null || !("kind" in p) || p.kind !== "job") return undefined;
    return "job_id" in p && typeof p.job_id === "string" ? p.job_id : undefined;
}

export function syncJobColumns(target: ColumnSyncTarget, server: ServerDataSource[]): void {
    for (const sds of server) {
        const ds = target.dataSources.find((d) => d.name === sds.name);
        if (!ds) continue;
        const store = ds.dataStore;
        const reload: string[] = [];
        for (const column of sds.columns ?? []) {
            const id = jobId(column);
            if (!id || column.deleted) continue;
            const local = store.columnIndex[column.field];
            if (!local) {
                store.addColumn(column);
            } else if (jobId(local) !== id && store.replaceColumnMetadata(column)) {
                reload.push(column.field);
            }
        }
        if (reload.length) {
            // not dirty: the new values came from the backend
            target.loadColumnSet(reload, ds.name, () => store.dataChanged(reload, false));
        }
    }
}
