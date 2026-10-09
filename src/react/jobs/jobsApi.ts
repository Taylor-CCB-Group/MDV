/**
 * Client for the per-project jobs routes (ADR-0012 "HTTP surface").
 *
 * `root` is the project root from `useProject()`: "" in single-project mode, "/project/<id>" in
 * multi-project mode.
 */

/** How the selector renders a param; mirrors `ParamSpec.type` in `mdvtools/jobs/registry.py`. */
export type JobParamType = "datasource" | "subgroup" | "column" | "dropdown" | "text" | "int" | "float";

export type JobParamValue = string | number;

export type JobParamSpec = {
    name: string;
    type: JobParamType;
    label: string;
    /** for "column" and "subgroup": the name of the datasource param whose columns or matrices are offered */
    options_from: string | null;
    default: JobParamValue | null;
    applies_to: string | null;
};

export type JobTool = {
    id: string;
    name: string;
    description: string;
    params: JobParamSpec[];
    output: { shape: string; datasource_param: string; columns_param: string };
    input_shape: string;
};

export type JobStatus =
    | "queued"
    | "staging"
    | "running"
    | "ingesting"
    | "done"
    | "failed"
    | "cancelled"
    | "stale"
    | "lost";

/** The client-safe view of a job record (`_client_view` in `server.py`). */
export type JobRecord = {
    job_id: string;
    tool_id: string;
    params: Record<string, JobParamValue>;
    status: JobStatus;
    input_filter_hash: string | null;
    /** seconds since the epoch */
    created: number;
    provenance: Record<string, unknown> | null;
    error: string | null;
};

/** Parse a JSON response, throwing the server's `{error}` message when the request failed. */
async function readJson(resp: Response) {
    if (resp.ok) return resp.json();
    // a Flask 500 is an HTML page, so the body may not be JSON
    const body = await resp.json().catch(() => null);
    throw new Error(typeof body?.error === "string" ? body.error : `${resp.status} ${resp.statusText}`);
}

export async function fetchTools(root: string): Promise<JobTool[]> {
    return readJson(await fetch(`${root}/jobs/tools`));
}

export async function submitJob(root: string, toolId: string, params: Record<string, JobParamValue>): Promise<string> {
    const resp = await fetch(`${root}/jobs`, {
        method: "POST",
        body: JSON.stringify({ tool_id: toolId, params }),
        headers: { "Content-Type": "application/json" },
    });
    const body: { job_id: string } = await readJson(resp);
    return body.job_id;
}

/** Every job of the project, newest first. */
export async function fetchJobs(root: string): Promise<JobRecord[]> {
    // GET /jobs returns records in no fixed order
    const jobs: JobRecord[] = await readJson(await fetch(`${root}/jobs`));
    return jobs.sort((a, b) => b.created - a.created);
}
