import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { renderHook, waitFor } from "@testing-library/react";
import type { ReactNode } from "react";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import type { ColumnSyncTarget, ServerDataSource } from "./jobColumns";
import type { JobRecord, JobStatus } from "./jobsApi";
import { useJobColumnSync } from "./useJobColumnSync";

const { syncJobColumns } = vi.hoisted(() => ({ syncJobColumns: vi.fn() }));
vi.mock("./jobColumns", () => ({ syncJobColumns }));

function job(job_id: string, status: JobStatus): JobRecord {
    return {
        job_id,
        tool_id: "umap",
        params: {},
        status,
        input_filter_hash: null,
        created: 1_700_000_000,
        provenance: null,
        error: null,
    };
}

const DATASOURCES: ServerDataSource[] = [{ name: "cells", columns: [{ field: "UMAP_1" }] }];
const TARGET: ColumnSyncTarget = { dataSources: [], loadColumnSet: vi.fn() };

let jobs: JobRecord[];
let configsOk: boolean;
let fetchMock: ReturnType<typeof vi.fn>;

function json(body: unknown, status = 200) {
    return new Response(JSON.stringify(body), { status, headers: { "Content-Type": "application/json" } });
}

beforeEach(() => {
    jobs = [];
    configsOk = true;
    syncJobColumns.mockReset();
    fetchMock = vi.fn(async (url: string) => {
        if (url === "/jobs") return json(jobs);
        if (url === "/get_configs")
            return configsOk ? json({ datasources: DATASOURCES, state: {} }) : json({ error: "down" }, 500);
        return json({ error: "not found" }, 404);
    });
    vi.stubGlobal("fetch", fetchMock);
});

afterEach(() => {
    vi.unstubAllGlobals();
});

function renderSync() {
    const client = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    const wrapper = ({ children }: { children: ReactNode }) => (
        <QueryClientProvider client={client}>{children}</QueryClientProvider>
    );
    renderHook(() => useJobColumnSync("", TARGET), { wrapper });
    // stands in for the next poll
    return { poll: () => client.invalidateQueries({ queryKey: ["jobs", ""] }) };
}

function configFetches() {
    return fetchMock.mock.calls.filter(([url]) => url === "/get_configs").length;
}

describe("useJobColumnSync", () => {
    test("syncs the page's columns when a job finishes", async () => {
        jobs = [job("a", "running")];
        const { poll } = renderSync();
        await waitFor(() => expect(fetchMock).toHaveBeenCalledWith("/jobs"));
        expect(syncJobColumns).not.toHaveBeenCalled();

        jobs = [job("a", "done")];
        await poll();

        await waitFor(() => expect(syncJobColumns).toHaveBeenCalledWith(TARGET, DATASOURCES));
    });

    test("syncs once per finished job", async () => {
        jobs = [job("a", "done")];
        const { poll } = renderSync();
        await waitFor(() => expect(syncJobColumns).toHaveBeenCalledTimes(1));

        jobs = [job("a", "done"), job("b", "running")];
        await poll();
        await waitFor(() => expect(fetchMock.mock.calls.filter(([url]) => url === "/jobs")).toHaveLength(2));

        expect(configFetches()).toBe(1);
        expect(syncJobColumns).toHaveBeenCalledTimes(1);
    });

    test("tries a failed sync again the next time the jobs are fetched", async () => {
        configsOk = false;
        jobs = [job("a", "done")];
        const { poll } = renderSync();
        await waitFor(() => expect(configFetches()).toBe(1));
        expect(syncJobColumns).not.toHaveBeenCalled();

        configsOk = true;
        await poll();

        await waitFor(() => expect(syncJobColumns).toHaveBeenCalledWith(TARGET, DATASOURCES));
    });
});
