import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import JobsPanel, { jobsRefetchInterval } from "./JobsPanel";
import type { JobsDataSource } from "./jobForm";
import type { JobRecord, JobTool } from "./jobsApi";

const TOOLS: JobTool[] = [
    {
        id: "umap",
        name: "UMAP",
        description: "Compute a UMAP embedding",
        input_shape: "matrix",
        output: { shape: "column", datasource_param: "datasource", columns_param: "output_name" },
        params: [
            {
                name: "datasource",
                type: "datasource",
                label: "Datasource",
                options_from: null,
                default: null,
                applies_to: null,
            },
            {
                name: "layer",
                type: "subgroup",
                label: "Matrix",
                options_from: "datasource",
                default: "gs",
                applies_to: null,
            },
            {
                name: "output_name",
                type: "text",
                label: "New column base name",
                options_from: null,
                default: "UMAP",
                applies_to: null,
            },
            {
                name: "n_neighbors",
                type: "int",
                label: "Neighbors",
                options_from: null,
                default: 15,
                applies_to: "neighbors",
            },
        ],
    },
    {
        id: "concat_columns",
        name: "Concatenate Columns",
        description: "Concatenate two columns",
        input_shape: "columns",
        output: { shape: "column", datasource_param: "datasource", columns_param: "output_name" },
        params: [
            {
                name: "datasource",
                type: "datasource",
                label: "Datasource",
                options_from: null,
                default: null,
                applies_to: null,
            },
            {
                name: "column_a",
                type: "column",
                label: "First column",
                options_from: "datasource",
                default: null,
                applies_to: null,
            },
            {
                name: "output_name",
                type: "text",
                label: "New column name",
                options_from: null,
                default: null,
                applies_to: null,
            },
        ],
    },
];

const SOURCES: JobsDataSource[] = [
    { name: "cells", columns: [{ field: "sample", name: "Sample" }], subgroups: ["gs"] },
];

function job(p: Partial<JobRecord> & Pick<JobRecord, "job_id">): JobRecord {
    return {
        tool_id: "umap",
        params: {},
        status: "done",
        input_filter_hash: null,
        created: 1_700_000_000,
        provenance: null,
        error: null,
        ...p,
    };
}

type Server = { jobs: JobRecord[]; submit: { status: number; body: unknown } };
let server: Server;
let fetchMock: ReturnType<typeof vi.fn>;

function json(body: unknown, status = 200) {
    return new Response(JSON.stringify(body), { status, headers: { "Content-Type": "application/json" } });
}

beforeEach(() => {
    server = { jobs: [], submit: { status: 202, body: { job_id: "new" } } };
    fetchMock = vi.fn(async (url: string, init?: RequestInit) => {
        if (url === "/jobs/tools") return json(TOOLS);
        if (url === "/jobs" && init?.method === "POST") {
            if (server.submit.status === 202)
                server.jobs.push(job({ job_id: "new", status: "queued", created: 2_000_000_000 }));
            return json(server.submit.body, server.submit.status);
        }
        if (url === "/jobs") return json(server.jobs);
        return json({ error: "not found" }, 404);
    });
    vi.stubGlobal("fetch", fetchMock);
});

afterEach(() => {
    vi.unstubAllGlobals();
});

function renderPanel() {
    const client = new QueryClient({ defaultOptions: { queries: { retry: false }, mutations: { retry: false } } });
    return render(
        <QueryClientProvider client={client}>
            <JobsPanel root="" sources={SOURCES} />
        </QueryClientProvider>,
    );
}

function inputValue(el: HTMLElement) {
    return el instanceof HTMLInputElement || el instanceof HTMLSelectElement ? el.value : undefined;
}

function postedBody() {
    const call = fetchMock.mock.calls.find(([, init]) => init?.method === "POST");
    return call ? JSON.parse(call[1].body) : undefined;
}

describe("JobsPanel", () => {
    test("shows the first tool's fields with their starting values", async () => {
        renderPanel();

        expect(inputValue(await screen.findByLabelText("Tool"))).toBe("umap");
        expect(screen.getByText("Compute a UMAP embedding")).toBeTruthy();
        expect(inputValue(screen.getByLabelText("Datasource"))).toBe("cells");
        expect(inputValue(screen.getByLabelText("Matrix"))).toBe("gs");
        expect(inputValue(screen.getByLabelText("New column base name"))).toBe("UMAP");
        expect(inputValue(screen.getByLabelText("Neighbors"))).toBe("15");
    });

    test("choosing another tool shows that tool's fields", async () => {
        renderPanel();

        fireEvent.change(await screen.findByLabelText("Tool"), { target: { value: "concat_columns" } });

        expect(inputValue(screen.getByLabelText("First column"))).toBe("");
        expect(screen.queryByLabelText("Matrix")).toBeNull();
    });

    test("running a job posts the form values and lists the new job", async () => {
        renderPanel();
        fireEvent.change(await screen.findByLabelText("Neighbors"), { target: { value: "30" } });

        fireEvent.click(screen.getByRole("button", { name: "Run" }));

        await waitFor(() =>
            expect(postedBody()).toEqual({
                tool_id: "umap",
                params: { datasource: "cells", layer: "gs", output_name: "UMAP", n_neighbors: 30 },
            }),
        );
        const row = await screen.findByRole("row", { name: /queued/ });
        expect(within(row).getByText("UMAP")).toBeTruthy();
    });

    test("a rejected submit shows the server's message", async () => {
        server.submit = { status: 400, body: { error: "'cells' is not a datasource" } };
        renderPanel();
        await screen.findByLabelText("Tool");

        fireEvent.click(screen.getByRole("button", { name: "Run" }));

        expect((await screen.findByRole("alert")).textContent).toContain("'cells' is not a datasource");
    });

    test("a numeric param that does not parse shows a message and sends nothing", async () => {
        renderPanel();
        fireEvent.change(await screen.findByLabelText("Neighbors"), { target: { value: "2.5" } });

        fireEvent.click(screen.getByRole("button", { name: "Run" }));

        expect((await screen.findByRole("alert")).textContent).toContain("Neighbors must be a whole number");
        expect(postedBody()).toBeUndefined();
    });

    test("lists jobs newest first, with the error of a failed job", async () => {
        server.jobs = [
            job({ job_id: "a", tool_id: "concat_columns", status: "done", created: 1_000 }),
            job({ job_id: "b", tool_id: "umap", status: "failed", created: 2_000, error: "worker exited 1" }),
        ];
        renderPanel();

        await screen.findByText("worker exited 1");
        const rows = screen.getAllByRole("row").slice(1); // skip the header row
        expect(within(rows[0]).getByText("failed")).toBeTruthy();
        expect(within(rows[1]).getByText("Concatenate Columns")).toBeTruthy();
    });
});

describe("jobsRefetchInterval", () => {
    test("polls while any job is still moving", () => {
        expect(jobsRefetchInterval([job({ job_id: "a", status: "running" })])).toBeGreaterThan(0);
        expect(jobsRefetchInterval([job({ job_id: "a", status: "queued" })])).toBeGreaterThan(0);
    });

    test("stops polling when every job has finished", () => {
        expect(
            jobsRefetchInterval([job({ job_id: "a", status: "done" }), job({ job_id: "b", status: "failed" })]),
        ).toBe(false);
        expect(jobsRefetchInterval(undefined)).toBe(false);
    });
});
