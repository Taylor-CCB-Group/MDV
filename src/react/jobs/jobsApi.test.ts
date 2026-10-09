import { afterEach, describe, expect, test, vi } from "vitest";
import { fetchJobs, fetchTools, submitJob } from "./jobsApi";

function jsonResponse(body: unknown, status = 200) {
    return new Response(JSON.stringify(body), {
        status,
        headers: { "Content-Type": "application/json" },
    });
}

function stubFetch(body: unknown, status = 200) {
    const fetchMock = vi.fn(async () => jsonResponse(body, status));
    vi.stubGlobal("fetch", fetchMock);
    return fetchMock;
}

afterEach(() => {
    vi.unstubAllGlobals();
});

describe("jobsApi", () => {
    test("fetchTools reads the tool registry from the project root", async () => {
        const tools = [{ id: "umap", name: "UMAP", params: [] }];
        const fetchMock = stubFetch(tools);

        await expect(fetchTools("/project/7")).resolves.toEqual(tools);
        expect(fetchMock).toHaveBeenCalledWith("/project/7/jobs/tools");
    });

    test("submitJob posts tool_id and params and returns the job id", async () => {
        const fetchMock = stubFetch({ job_id: "abc" }, 202);

        await expect(submitJob("", "umap", { datasource: "cells" })).resolves.toBe("abc");

        expect(fetchMock).toHaveBeenCalledWith(
            "/jobs",
            expect.objectContaining({
                method: "POST",
                body: JSON.stringify({ tool_id: "umap", params: { datasource: "cells" } }),
            }),
        );
    });

    test("submitJob throws the server's error message on a rejected submit", async () => {
        stubFetch({ error: "'nope' is not a datasource" }, 400);

        await expect(submitJob("", "umap", { datasource: "nope" })).rejects.toThrow("'nope' is not a datasource");
    });

    test("a non-JSON error response throws the HTTP status", async () => {
        // a Flask 500 is an HTML page, not {error}
        vi.stubGlobal(
            "fetch",
            vi.fn(async () => new Response("<html>boom</html>", { status: 500 })),
        );

        await expect(fetchTools("")).rejects.toThrow("500");
    });

    test("fetchJobs returns records newest first", async () => {
        // GET /jobs returns records in no fixed order
        stubFetch([
            { job_id: "old", created: 100 },
            { job_id: "new", created: 300 },
            { job_id: "mid", created: 200 },
        ]);

        const jobs = await fetchJobs("");

        expect(jobs.map((j) => j.job_id)).toEqual(["new", "mid", "old"]);
    });
});
