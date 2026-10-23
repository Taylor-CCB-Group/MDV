import { useQuery } from "@tanstack/react-query";
import { useEffect, useRef } from "react";
import { type ColumnSyncTarget, syncJobColumns } from "./jobColumns";
import { type JobRecord, type JobStatus, fetchConfigs, fetchJobs } from "./jobsApi";

/** Statuses the driver will still move on from (ADR-0012), so the list keeps polling. */
const MOVING: JobStatus[] = ["queued", "staging", "running", "ingesting"];
const POLL_MS = 2000;

export function jobsRefetchInterval(jobs: JobRecord[] | undefined): number | false {
    return jobs?.some((j) => MOVING.includes(j.status)) ? POLL_MS : false;
}

/** The project's jobs, polled while any is still moving; shared by every component that lists them. */
export function jobsQueryOptions(root: string) {
    return {
        queryKey: ["jobs", root],
        queryFn: () => fetchJobs(root),
        refetchInterval: (query: { state: { data: JobRecord[] | undefined } }) => jobsRefetchInterval(query.state.data),
    };
}

/**
 * Follows the project's jobs and, the first time each one is seen done, brings the columns it
 * wrote into the open page (`syncJobColumns`). Jobs already done when the page opens sync too;
 * the sync changes nothing for columns the page loaded with.
 */
export function useJobColumnSync(root: string, target: ColumnSyncTarget): void {
    // dataUpdatedAt changes on every fetch, even when the jobs did not, so a failed sync is retried
    const { data: jobs, dataUpdatedAt } = useQuery(jobsQueryOptions(root));
    const synced = useRef(new Set<string>());

    useEffect(() => {
        const fresh = (jobs ?? []).filter((j) => j.status === "done" && !synced.current.has(j.job_id));
        if (fresh.length === 0) return;
        fetchConfigs(root)
            .then((configs) => {
                syncJobColumns(target, configs.datasources);
                // only after a successful sync, so a failed one is tried again on the next fetch
                for (const j of fresh) synced.current.add(j.job_id);
            })
            .catch((e) => console.error("could not add the new job columns to the page", e));
    }, [jobs, dataUpdatedAt, root, target]);
}
